-- =============================================================================
-- holt_winters_q2_30_opt_v11.vhd
--
-- Builds on v10. All v10 changes retained.
--
-- v11 CHANGE -- THREE REDUNDANT STATES MERGED:
--
--   MERGE 1: S_VALIDATE eliminated.
--     Both checks (M<2, N<2*M) moved directly into the two S_COLLECT exit
--     branches (sample_cnt=MAX_N-1 and sample_cnt>0). On pass, S_COLLECT
--     jumps straight to S_LOAD_CONSTANTS. On fail, jumps to S_ERROR.
--     Saves 1 cycle per run.
--
--   MERGE 2: S_INIT_LT_DIV eliminated.
--     The divide (fp_mul_data on acc_L and acc_T) moved into the else branch
--     of S_INIT_LT, which was previously a dead redirect cycle. L_reg, T_reg,
--     and all downstream resets now fire in that same else. Saves 1 cycle.
--
--   MERGE 3: S_SAVG_DIV eliminated.
--     The mean divide moved into the else branch of S_SAVG, same pattern.
--     cur_season_mean, acc_L reset, flat_idx reset all fire there.
--     Saves MAX_SEASONS cycles (5 cycles in the 60-sample/12-season config).
--
-- Total states removed: 3  (S_VALIDATE, S_INIT_LT_DIV, S_SAVG_DIV)
-- Total cycles saved:   7  (1 + 1 + 5) for the default config.
-- FSM state count:      13 (was 16 in v10).
--
-- v11 ALSO -- FIRST-SAMPLE MECHANISM SIMPLIFIED:
--
--   REMOVED: first_sample_reg (32-bit FF register)
--   REMOVED: first_sample_vld (1-bit flag FF)
--   SAVING:  33 FFs
--
--   OLD MECHANISM (v9/v10):
--     S_IDLE captured data_in into first_sample_reg and set first_sample_vld.
--     S_COLLECT checked the flag on every cycle to decide between writing
--     first_sample_reg or data_in. Two separate branches, one flag register.
--
--   NEW MECHANISM (v11):
--     S_IDLE transitions to S_COLLECT with sample_cnt=0. No data captured.
--     S_COLLECT uses sample_cnt=0 as the implicit first-sample signal:
--       when sample_cnt=0 and valid_in='1': write data_in into data_buf(0).
--       This works because the top module (hw_top ROM_PREFETCH state) holds
--       data_in and valid_in stable for at least one cycle after S_IDLE fires.
--       So when S_COLLECT first executes, data_in still holds sample[0].
--     All subsequent cycles: normal stream, sample_cnt increments each write.
--
--   SINGLE WRITE CONDITION PRESERVED:
--     data_buf is ONLY ever written inside S_COLLECT via the single
--     `data_buf(sample_cnt) <= data_in` path. S_IDLE no longer touches
--     data_buf at all. Vivado sees one write port, one write enable,
--     one address bus. Distributed RAM inference is maintained. 5 replicas.
--
--   REQUIREMENT ON TOP MODULE:
--     The top module MUST hold data_in and valid_in stable for at least
--     2 rising edges after the first valid_in pulse. hw_top already does
--     this via ROM_PREFETCH (data pre-loaded before valid asserted) and
--     ROM_STREAM (valid held continuous). No top-module change needed.
--
-- All timing-critical paths (DSP pipeline, RAM prefetch) are unchanged.
-- WNS, resource counts, and distributed RAM inference are unaffected.
--
-- PORT SIMPLIFICATION (carried from v9):
--   start, forecast_start, m_in, horizon_in all removed.
--   M_act       <= MAX_M       (register, will be driven via AXI from PS later)
--   HORIZON_act <= MAX_HORIZON (register, will be driven via AXI from PS later)
--   SEASONS_act <= MAX_SEASONS (register, PS computes N/M and passes as generic)
--   FSM fully automatic -- driven only by valid_in stream.
-- =============================================================================

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

-- =============================================================================
-- Generics Declaration -- Values from top_module are used
-- =============================================================================
entity holt_winters_q2_30 is
    generic (
        MAX_N       : integer := 72;  -- Maximum sample buffer depth; sets data_buf array size and all loop upper bounds
        MAX_M       : integer := 24;  -- Maximum season length; sets S_arr, S_acc array sizes and ROM upper index
        MAX_HORIZON : integer := 24;  -- Maximum forecast steps; controls S_FORECAST loop length
        MAX_SEASONS : integer := 8    -- PS sets this to actual N/M (e.g. 60/12=5); bounds SEASONS_act and S_acc accumulation
    );
    -- ==========================================================================
    -- Ports Declared
    -- ==========================================================================
    port (
        clk            : in  std_logic;                    -- Standard synchronous clock
        rst            : in  std_logic;                    -- Active high reset
        alpha          : in  signed(31 downto 0);          -- Q2.30 smoothing parameter for level
        beta           : in  signed(31 downto 0);          -- Q2.30 smoothing parameter for trend
        gamma          : in  signed(31 downto 0);          -- Q2.30 smoothing parameter for seasonality
        data_in        : in  signed(31 downto 0);          -- Q16.16 input sample
        valid_in       : in  std_logic;                    -- Handshake signal, ALSO the FSM trigger replacing start/forecast_start
        fitted_out     : out signed(31 downto 0);          -- Reconstructed training values
        fitted_valid   : out std_logic;                    -- Strobes high for one cycle when fitted_out is valid
        forecast_out   : out signed(31 downto 0);          -- Prediction output
        forecast_valid : out std_logic;                    -- Strobes high for one cycle when forecast_out is valid
        last_forecast  : out std_logic;                    -- Pulses high on the final forecast step, signals horizon reached
        error_out      : out std_logic;                    -- Pulses high for one cycle when a validation error is detected
        error_code     : out std_logic_vector(2 downto 0)  -- 3-bit code identifying which validation check failed
    );
end entity holt_winters_q2_30;

architecture rtl of holt_winters_q2_30 is

    subtype fp_data  is signed(31 downto 0); -- Q16.16: 16-bit integer, 16-bit fraction
    subtype fp_param is signed(31 downto 0); -- Q2.30: 2-bit integer, 30-bit fraction (for 0..1 params)

    -- 1.0 in Q2.30 = 2^30; subtracted from alpha/beta/gamma in S_UPDATE and S_UPDATE_WAIT
    -- to compute (1-alpha), (1-beta), (1-gamma) without a hardware subtractor on the parameter bus.
    constant ONE_PARAM : fp_param := to_signed(1073741824, 32);

    signal N_act       : integer range 1 to MAX_N       := 1; -- Actual number of samples received; set in S_COLLECT when stream ends
    signal M_act       : integer range 2 to MAX_M       := 2; -- Season length; assigned from MAX_M in S_IDLE, will be driven via AXI later
    signal HORIZON_act : integer range 1 to MAX_HORIZON := 1; -- Forecast steps; assigned from MAX_HORIZON in S_LOAD_CONSTANTS, will be AXI later

    -- SEASONS_act: register holding the actual season count.
    -- Currently loaded from MAX_SEASONS generic (PS computes N/M and passes it).
    -- Initialised to MAX_SEASONS so it is valid immediately on reset.
    -- Will be driven via AXI from PS later.
    signal SEASONS_act : integer range 1 to MAX_SEASONS := MAX_SEASONS;

    type data_array   is array (0 to MAX_N-1) of fp_data; -- Ring of MAX_N samples; distributed LUT RAM
    type season_array is array (0 to MAX_M-1) of fp_data; -- Holds M seasonal factors S[0..M-1]

    -- =========================================================================
    -- data_buf: distributed RAM. CRITICAL: written ONLY inside S_COLLECT
    -- with sample_cnt as the sole address. Any additional write condition
    -- outside S_COLLECT breaks Vivado RAM inference and causes 1920 FF
    -- fallback (+2000 LUTs, +1900 FFs vs baseline).
    -- S_IDLE does NOT write data_buf. sample_cnt=0 inside S_COLLECT is used
    -- as the implicit first-sample trigger instead of a separate flag register.
    -- =========================================================================
    signal data_buf : data_array := (others => (others => '0'));
    attribute ram_style : string;
    attribute ram_style of data_buf : signal is "distributed"; -- Forces distributed LUT RAM inference

    -- =========================================================================
    -- first_sample_reg and first_sample_vld REMOVED in v11.
    -- S_COLLECT now uses sample_cnt=0 as the implicit first-sample condition.
    -- The top module holds data_in and valid_in stable so sample[0] is still
    -- present on data_in when S_COLLECT executes its first cycle.
    -- =========================================================================

    signal S_arr           : season_array := (others => (others => '0')); -- Current seasonal factors S[0..M-1], updated each training cycle
    signal S_acc           : season_array := (others => (others => '0')); -- Accumulator for seasonal factor sums during initialisation
    signal cur_season_mean : fp_data      := (others => '0');             -- Mean of one season window; held for S_SACC subtraction

    signal L_reg : fp_data := (others => '0'); -- Current level estimate L_t; committed in S_UPDATE_WAIT
    signal T_reg : fp_data := (others => '0'); -- Current trend estimate T_t; committed in S_UPDATE_FINALIZE

    -- =========================================================================
    -- Reciprocal ROM: stores floor(65536/i) for i=1..MAX_M^2 in Q16.16 format.
    -- Replaces all hardware division with a single-cycle multiply.
    -- Built at elaboration time; synthesises to LUT ROM.
    -- =========================================================================
    type rcp_rom_t is array (1 to MAX_M * MAX_M) of fp_data;

    function init_rcp_rom return rcp_rom_t is
        variable rom : rcp_rom_t;
    begin
        for i in 1 to MAX_M * MAX_M loop
            rom(i) := resize(to_signed((65536 + i/2) / i, 32), 32); -- Rounded reciprocal in Q16.16
        end loop;
        return rom;
    end function;

    constant RCP_ROM : rcp_rom_t := init_rcp_rom;

    signal inv_M_reg  : fp_data := (others => '0'); -- 1/M in Q16.16; loaded in S_LOAD_CONSTANTS; used in S_INIT_LT and S_SAVG else branches
    signal inv_M2_reg : fp_data := (others => '0'); -- 1/M^2 in Q16.16; loaded in S_LOAD_CONSTANTS; used in S_INIT_LT else for trend init
    signal inv_S_reg  : fp_data := (others => '0'); -- 1/SEASONS in Q16.16; loaded in S_LOAD_CONSTANTS; used in S_SACC normalisation

    signal M_act_m1         : integer range 0 to MAX_M-1       := 1; -- M-1; precomputed to avoid a subtractor on si_counter wrap check every cycle
    signal SEASONS_m1       : integer range 0 to MAX_SEASONS-1 := 0; -- SEASONS-1; precomputed for s_cnt loop termination in S_SACC
    signal next_season_flat : integer range 0 to MAX_N-1       := 0; -- Precomputed season_base + M_act; avoids adder in S_SACC hot path

    -- =========================================================================
    -- FSM state declarations.
    -- v11: S_VALIDATE, S_INIT_LT_DIV, S_SAVG_DIV removed (merged into neighbours).
    -- =========================================================================
    type state_t is (
        S_IDLE,            -- Waits for first valid_in; sets M_act and sample_cnt=0; transitions to S_COLLECT
        S_COLLECT,         -- Sole writer of data_buf; sample_cnt=0 captures sample[0]; inline validation on exit
        S_LOAD_CONSTANTS,  -- One-cycle ROM prefetch; loads inv_M_reg, inv_M2_reg, inv_S_reg; precomputes helper signals
        S_INIT_LT,         -- Accumulates level and trend sums over first season; divides in else branch (merged S_INIT_LT_DIV)
        S_SAVG,            -- Accumulates one season window into acc_L; divides in else branch (merged S_SAVG_DIV)
        S_SACC,            -- Subtracts cur_season_mean from each element; accumulates into S_acc; normalises when all seasons done
        S_UPDATE_PRE,      -- Prefetches S_arr(si_counter) and data_buf(t_cnt) one cycle early to hide distributed RAM read latency
        S_UPDATE,          -- Computes fitted value and launches alpha DSP multiplies
        S_UPDATE_MULT,     -- Collects alpha DSP results to form L_new; re-latches data sample for gamma term
        S_UPDATE_WAIT,     -- Emits fitted_out; commits L_reg; launches beta and gamma DSP multiplies
        S_UPDATE_FINALIZE, -- Commits T_reg and S_arr update; loops back to S_UPDATE or advances to S_FORECAST
        S_FORECAST,        -- Emits one forecast per cycle for HORIZON_act cycles then returns to S_IDLE
        S_ERROR            -- Pulses error_out for one cycle then returns to S_IDLE
    );
    signal state : state_t := S_IDLE;

    -- =========================================================================
    -- Counters
    -- =========================================================================
    signal sample_cnt : integer range 0 to MAX_N       := 0; -- Write pointer in S_COLLECT; 0 on first entry = implicit first-sample trigger
    signal i_cnt      : integer range 0 to MAX_M       := 0; -- Inner loop index in S_INIT_LT, S_SAVG, S_SACC
    signal s_cnt      : integer range 0 to MAX_SEASONS := 0; -- Season loop index in S_SACC
    signal t_cnt      : integer range 0 to MAX_N       := 0; -- Training time index in UPDATE pipeline; incremented in S_UPDATE
    signal k_cnt      : integer range 0 to MAX_HORIZON := 0; -- Forecast step counter; steps 0 to HORIZON_act-1 in S_FORECAST

    signal flat_idx    : integer range 0 to MAX_N-1 := 0; -- Linear data_buf address in S_SAVG and S_SACC
    signal season_base : integer range 0 to MAX_N-1 := 0; -- Start address of current season window in data_buf

    signal acc_L : fp_data := (others => '0'); -- Running sum for level init (S_INIT_LT) and season mean (S_SAVG)
    signal acc_T : fp_data := (others => '0'); -- Running sum of season-over-season differences in S_INIT_LT

    -- =========================================================================
    -- UPDATE pipeline registers
    -- =========================================================================
    signal upd_S_pre     : fp_data := (others => '0'); -- S[si] latched one cycle early in S_UPDATE_PRE to hide RAM latency
    signal upd_data_pre  : fp_data := (others => '0'); -- data_buf[t] latched one cycle early in S_UPDATE_PRE
    signal upd_data_wait : fp_data := (others => '0'); -- data_buf[upd_t] re-latched in S_UPDATE_MULT for gamma term
    signal prev_L        : fp_data := (others => '0'); -- L_reg before overwrite; needed for beta term (L_new - L_old)

    signal upd_si    : integer range 0 to MAX_M-1 := 0; -- Season index forwarded so S_UPDATE_FINALIZE writes S_arr correctly
    signal upd_fcast : fp_data   := (others => '0');    -- Fitted value (L+T+S) held from S_UPDATE until S_UPDATE_WAIT
    signal upd_S_cur : fp_data   := (others => '0');    -- S[si] copy forwarded for gamma term in S_UPDATE_WAIT
    signal upd_newL  : fp_data   := (others => '0');    -- L_new formed in S_UPDATE_MULT; consumed in S_UPDATE_WAIT

    -- Frozen sample index; t_cnt increments in S_UPDATE so upd_t preserves the
    -- current index for the deferred re-read in S_UPDATE_MULT and the termination
    -- check in S_UPDATE_FINALIZE.
    signal upd_t : integer range 0 to MAX_N-1 := 0;

    signal pipe_valid : std_logic := '0'; -- High when pipe_newL is valid; S_UPDATE uses pipe_newL+T_reg instead of stale L_reg+T_reg
    signal pipe_newL  : fp_data   := (others => '0'); -- Forwards L_new before L_reg is committed; prevents stale-level error

    -- =========================================================================
    -- DSP pipeline registers -- one per multiply term in the three HW equations
    -- =========================================================================
    signal mul_alpha_tmp : fp_data := (others => '0'); -- alpha    * (x - S)
    signal mul_1ma_tmp   : fp_data := (others => '0'); -- (1-alpha)* (L + T)
    signal mul_beta_tmp  : fp_data := (others => '0'); -- beta     * (L_new - L_old)
    signal mul_1mb_tmp   : fp_data := (others => '0'); -- (1-beta) * T
    signal mul_gamma_tmp : fp_data := (others => '0'); -- gamma    * (x - L_new)
    signal mul_1mg_tmp   : fp_data := (others => '0'); -- (1-gamma)* S

    signal si_counter    : integer range 0 to MAX_M-1 := 0; -- Season index for current training sample; wraps modulo M
    signal fc_si_counter : integer range 0 to MAX_M-1 := 0; -- Season index for forecast; starts one ahead of si_counter

    -- =========================================================================
    -- fp_mul_data: Q16.16 x Q16.16 -> Q16.16
    -- =========================================================================
    function fp_mul_data(a, b : fp_data) return fp_data is
        variable tmp : signed(63 downto 0);
        attribute use_dsp : string;
        attribute use_dsp of tmp : variable is "yes";
    begin
        tmp := a * b + resize(to_signed(32768, 32), 64); -- Round to nearest at Q16.16 boundary
        return tmp(47 downto 16);
    end function;

    -- =========================================================================
    -- fp_mul_param: Q2.30 x Q16.16 -> Q16.16
    -- =========================================================================
    function fp_mul_param(param : fp_param; data : fp_data) return fp_data is
        variable tmp : signed(63 downto 0);
        attribute use_dsp : string;
        attribute use_dsp of tmp : variable is "yes";
    begin
        tmp := param * data + resize(to_signed(536870912, 32), 64); -- Round at Q2.30 boundary
        return tmp(61 downto 30);
    end function;

begin

    process(clk)
        variable v_1ma    : fp_param;              -- (1-alpha): computed in S_UPDATE
        variable v_1mb    : fp_param;              -- (1-beta):  computed in S_UPDATE_WAIT
        variable v_1mg    : fp_param;              -- (1-gamma): computed in S_UPDATE_WAIT
        variable k1_fp    : fp_data;               -- (k+1) in Q16.16 for trend projection in S_FORECAST
        variable tmp_data : fp_data;               -- (x - S): de-seasonalised input in S_UPDATE
        variable tmp_LpT  : fp_data;               -- (L + T): one-step-ahead prediction in S_UPDATE
        variable next_si  : integer range 0 to MAX_M-1; -- Next si_counter; computed in S_UPDATE_FINALIZE for immediate prefetch
    begin
    if rising_edge(clk) then

        -- Default all strobes low every cycle
        forecast_valid <= '0';
        fitted_valid   <= '0';
        last_forecast  <= '0';
        error_out      <= '0';
        error_code     <= (others => '0');

        if rst = '1' then
            state            <= S_IDLE;
            SEASONS_act      <= MAX_SEASONS;
            sample_cnt       <= 0;
            i_cnt            <= 0;
            s_cnt            <= 0;
            t_cnt            <= 0;
            k_cnt            <= 0;
            flat_idx         <= 0;
            season_base      <= 0;
            acc_L            <= (others => '0');
            acc_T            <= (others => '0');
            L_reg            <= (others => '0');
            T_reg            <= (others => '0');
            prev_L           <= (others => '0');
            pipe_valid       <= '0';
            pipe_newL        <= (others => '0');
            S_acc            <= (others => (others => '0'));
            S_arr            <= (others => (others => '0'));
            cur_season_mean  <= (others => '0');
            si_counter       <= 0;
            fc_si_counter    <= 0;
            upd_S_pre        <= (others => '0');
            upd_data_pre     <= (others => '0');
            upd_data_wait    <= (others => '0');
            upd_si           <= 0;
            upd_fcast        <= (others => '0');
            upd_S_cur        <= (others => '0');
            upd_newL         <= (others => '0');
            upd_t            <= 0;
            inv_M_reg        <= (others => '0');
            inv_M2_reg       <= (others => '0');
            inv_S_reg        <= (others => '0');
            M_act_m1         <= 1;
            SEASONS_m1       <= 0;
            next_season_flat <= 0;
            mul_alpha_tmp    <= (others => '0');
            mul_1ma_tmp      <= (others => '0');
            mul_beta_tmp     <= (others => '0');
            mul_1mb_tmp      <= (others => '0');
            mul_gamma_tmp    <= (others => '0');
            mul_1mg_tmp      <= (others => '0');
            -- Note: first_sample_reg and first_sample_vld removed in v11

        else
            case state is

                -- =============================================================
                -- S_IDLE: waits for first valid_in pulse.
                --
                -- On valid_in='1': set M_act, reset sample_cnt to 0,
                -- transition to S_COLLECT. No data captured here.
                --
                -- data_buf is NOT written in S_IDLE. Single write condition
                -- preserved. sample_cnt=0 inside S_COLLECT acts as the
                -- implicit trigger to capture sample[0] from data_in.
                --
                -- REQUIREMENT: top module must hold data_in and valid_in
                -- stable for at least the first S_COLLECT cycle so that
                -- data_in still presents sample[0] when S_COLLECT executes.
                -- hw_top ROM_PREFETCH + ROM_STREAM already guarantees this.
                -- =============================================================
                when S_IDLE =>
                    if valid_in = '1' then
                        M_act      <= MAX_M; -- Load season length; will come from AXI later
                        sample_cnt <= 0;     -- Reset write pointer; sample_cnt=0 in S_COLLECT captures sample[0]
                        state      <= S_COLLECT;
                    end if;

                -- =============================================================
                -- S_COLLECT: sole and only writer of data_buf.
                --
                -- FIRST-SAMPLE MECHANISM (v11, no flag register):
                --   On first entry sample_cnt=0 and valid_in='1' (top holds stream).
                --   data_in still holds sample[0] because top module keeps it stable.
                --   The normal `valid_in='1'` branch writes data_in to data_buf(0).
                --   sample_cnt increments to 1. Next cycle catches sample[1] normally.
                --   No separate branch, no flag signal -- sample_cnt=0 is the implicit key.
                --
                -- WRITE CONDITION:
                --   data_buf written only here via data_buf(sample_cnt) <= data_in.
                --   Vivado sees one write port, one address, one write enable.
                --   Distributed RAM inference maintained. 5 replicas preserved.
                --
                -- EXIT (inline validation merged from S_VALIDATE):
                --   Buffer full (sample_cnt = MAX_N-1): validate then S_LOAD_CONSTANTS.
                --   Stream ended early (valid_in dropped): validate then S_LOAD_CONSTANTS.
                --   On fail: S_ERROR with error_code.
                -- =============================================================
                when S_COLLECT =>
                    if valid_in = '1' then
                        data_buf(sample_cnt) <= data_in; -- Write sample[sample_cnt]; when cnt=0 this writes sample[0]

                        if sample_cnt = MAX_N - 1 then
                            -- Buffer full; inline validation
                            if MAX_M < 2 then
                                error_out  <= '1';
                                error_code <= "010"; -- M too small
                                state      <= S_ERROR;
                            elsif MAX_N < 2 * MAX_M then
                                error_out  <= '1';
                                error_code <= "001"; -- Not enough samples for two seasons
                                state      <= S_ERROR;
                            else
                                N_act      <= MAX_N;           -- Lock in sample count
                                sample_cnt <= 0;               -- Reset pointer
                                acc_L      <= (others => '0'); -- Clear for S_INIT_LT
                                acc_T      <= (others => '0'); -- Clear for S_INIT_LT
                                i_cnt      <= 0;
                                state      <= S_LOAD_CONSTANTS;
                            end if;
                        else
                            sample_cnt <= sample_cnt + 1; -- Advance write pointer
                        end if;

                    elsif sample_cnt > 0 then
                        -- valid_in dropped before buffer full; sample_cnt = actual N
                        if MAX_M < 2 then
                            error_out  <= '1';
                            error_code <= "010";
                            state      <= S_ERROR;
                        elsif sample_cnt < 2 * MAX_M then
                            error_out  <= '1';
                            error_code <= "001";
                            state      <= S_ERROR;
                        else
                            N_act      <= sample_cnt;      -- Actual sample count
                            sample_cnt <= 0;
                            acc_L      <= (others => '0');
                            acc_T      <= (others => '0');
                            i_cnt      <= 0;
                            state      <= S_LOAD_CONSTANTS;
                        end if;
                    end if;

                -- =============================================================
                -- S_LOAD_CONSTANTS: one-cycle ROM prefetch.
                -- All reciprocals and precomputed counters loaded and registered
                -- here so they are stable when S_INIT_LT reads them next cycle.
                -- =============================================================
                when S_LOAD_CONSTANTS =>
                    HORIZON_act      <= MAX_HORIZON;              -- Load forecast horizon; will come from AXI later
                    inv_M_reg        <= RCP_ROM(M_act);           -- 1/M for level init and season mean
                    inv_M2_reg       <= RCP_ROM(M_act * M_act);   -- 1/M^2 for trend init
                    inv_S_reg        <= RCP_ROM(SEASONS_act);     -- 1/SEASONS for seasonal normalisation
                    M_act_m1         <= M_act - 1;                -- M-1 for si_counter wrap
                    SEASONS_m1       <= SEASONS_act - 1;          -- SEASONS-1 for s_cnt termination
                    next_season_flat <= M_act;                    -- First inter-season boundary
                    i_cnt            <= 0;
                    state            <= S_INIT_LT;

                -- =============================================================
                -- S_INIT_LT: accumulates level and trend sums over first season.
                -- acc_L += x[i]            -> L0 = mean of first season
                -- acc_T += x[M+i] - x[i]  -> T0 = mean trend per step
                -- Runs M cycles then divides in else branch (merged S_INIT_LT_DIV).
                -- =============================================================
                when S_INIT_LT =>
                    if i_cnt < M_act then
                        acc_L <= acc_L + data_buf(i_cnt);                              -- Accumulate x[i]
                        acc_T <= acc_T + (data_buf(M_act + i_cnt) - data_buf(i_cnt)); -- Accumulate x[M+i]-x[i]
                        i_cnt <= i_cnt + 1;
                    else
                        -- Divide (was S_INIT_LT_DIV)
                        L_reg       <= fp_mul_data(acc_L, inv_M_reg);  -- L0 = acc_L / M
                        T_reg       <= fp_mul_data(acc_T, inv_M2_reg); -- T0 = acc_T / M^2
                        S_acc       <= (others => (others => '0'));     -- Clear seasonal accumulator
                        acc_L       <= (others => '0');
                        acc_T       <= (others => '0');
                        i_cnt       <= 0;
                        s_cnt       <= 0;
                        flat_idx    <= 0;
                        season_base <= 0;
                        state       <= S_SAVG;
                    end if;

                -- =============================================================
                -- S_SAVG: accumulates one season window into acc_L.
                -- Runs M cycles stepping flat_idx through data_buf.
                -- Divides in else branch to get cur_season_mean (merged S_SAVG_DIV).
                -- =============================================================
                when S_SAVG =>
                    if i_cnt < M_act then
                        acc_L    <= acc_L + data_buf(flat_idx); -- Accumulate season element
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                    else
                        -- Divide (was S_SAVG_DIV)
                        cur_season_mean <= fp_mul_data(acc_L, inv_M_reg); -- mean = acc_L / M
                        acc_L           <= (others => '0');
                        i_cnt           <= 0;
                        flat_idx        <= season_base; -- Reset to season start for S_SACC
                        state           <= S_SACC;
                    end if;

                -- =============================================================
                -- S_SACC: subtracts cur_season_mean from each element;
                -- accumulates into S_acc. After all seasons normalises into S_arr.
                -- =============================================================
                when S_SACC =>
                    if i_cnt < M_act then
                        S_acc(i_cnt) <= S_acc(i_cnt) +
                                        data_buf(flat_idx) - cur_season_mean; -- (x - mean)
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                    else
                        if s_cnt < SEASONS_m1 then
                            -- More seasons; advance window and loop to S_SAVG
                            s_cnt            <= s_cnt + 1;
                            i_cnt            <= 0;
                            season_base      <= next_season_flat;
                            flat_idx         <= next_season_flat;
                            next_season_flat <= next_season_flat + M_act;
                            state            <= S_SAVG;
                        else
                            -- All seasons done; normalise into S_arr
                            for j in 0 to MAX_M - 1 loop
                                if j < M_act then
                                    S_arr(j) <= fp_mul_data(S_acc(j), inv_S_reg); -- S[j] / SEASONS
                                end if;
                            end loop;
                            si_counter    <= 0;
                            fc_si_counter <= 0;
                            i_cnt         <= 0;
                            s_cnt         <= 0;
                            t_cnt         <= 0;
                            pipe_valid    <= '0';
                            state         <= S_UPDATE_PRE;
                        end if;
                    end if;

                -- =============================================================
                -- S_UPDATE_PRE: one-cycle RAM prefetch before S_UPDATE.
                -- Latches S_arr(si_counter) and data_buf(t_cnt) one cycle early
                -- to hide distributed RAM read latency.
                -- =============================================================
                when S_UPDATE_PRE =>
                    upd_S_pre    <= S_arr(si_counter);  -- Seasonal factor for this time step
                    upd_data_pre <= data_buf(t_cnt);    -- Raw sample for this time step
                    upd_si       <= si_counter;          -- Freeze season index for pipeline
                    upd_t        <= t_cnt;               -- Freeze sample index for pipeline
                    state        <= S_UPDATE;

                -- =============================================================
                -- S_UPDATE: launches alpha DSP multiplies.
                -- =============================================================
                when S_UPDATE =>
                    v_1ma    := ONE_PARAM - alpha;         -- (1-alpha)
                    tmp_data := upd_data_pre - upd_S_pre;  -- x - S

                    if pipe_valid = '1' then tmp_LpT := pipe_newL + T_reg; -- Warm pipeline
                    else                     tmp_LpT := L_reg    + T_reg;  -- Cold start
                    end if;

                    upd_fcast     <= tmp_LpT + upd_S_pre;           -- Fitted = (L+T)+S
                    upd_S_cur     <= upd_S_pre;                      -- Forward S for gamma term
                    mul_alpha_tmp <= fp_mul_param(alpha, tmp_data);  -- DSP: alpha*(x-S)
                    mul_1ma_tmp   <= fp_mul_param(v_1ma, tmp_LpT);  -- DSP: (1-alpha)*(L+T)
                    pipe_valid    <= '1';
                    t_cnt         <= upd_t + 1;                      -- Advance training index
                    state         <= S_UPDATE_MULT;

                -- =============================================================
                -- S_UPDATE_MULT: collects alpha DSP results; forms L_new.
                -- Re-latches x[upd_t] for gamma term (t_cnt has advanced).
                -- =============================================================
                when S_UPDATE_MULT =>
                    upd_newL      <= mul_alpha_tmp + mul_1ma_tmp; -- L_new
                    pipe_newL     <= mul_alpha_tmp + mul_1ma_tmp; -- Forward L_new before L_reg committed
                    upd_data_wait <= data_buf(upd_t);             -- Re-latch x[t] at frozen index
                    state         <= S_UPDATE_WAIT;

                -- =============================================================
                -- S_UPDATE_WAIT: emits fitted_out; launches beta/gamma DSPs.
                -- =============================================================
                when S_UPDATE_WAIT =>
                    v_1mb := ONE_PARAM - beta;  -- (1-beta)
                    v_1mg := ONE_PARAM - gamma; -- (1-gamma)

                    fitted_out   <= upd_fcast; -- Emit reconstructed value
                    fitted_valid <= '1';

                    prev_L <= L_reg;    -- Snapshot before overwrite
                    L_reg  <= upd_newL; -- Commit new level

                    mul_beta_tmp  <= fp_mul_param(beta,  upd_newL - prev_L);        -- DSP: beta*(L_new-L_old)
                    mul_1mb_tmp   <= fp_mul_param(v_1mb, T_reg);                    -- DSP: (1-beta)*T
                    mul_gamma_tmp <= fp_mul_param(gamma, upd_data_wait - upd_newL); -- DSP: gamma*(x-L_new)
                    mul_1mg_tmp   <= fp_mul_param(v_1mg, upd_S_cur);                -- DSP: (1-gamma)*S

                    pipe_valid <= '0';
                    state      <= S_UPDATE_FINALIZE;

                -- =============================================================
                -- S_UPDATE_FINALIZE: commits T_reg and S_arr; loops or forecasts.
                -- Inline prefetch skips S_UPDATE_PRE on all but the first iteration.
                -- =============================================================
                when S_UPDATE_FINALIZE =>
                    T_reg         <= mul_beta_tmp + mul_1mb_tmp;  -- New trend
                    S_arr(upd_si) <= mul_gamma_tmp + mul_1mg_tmp; -- Updated seasonal factor

                    if si_counter = M_act_m1 then next_si := 0;
                    else                          next_si := si_counter + 1;
                    end if;
                    si_counter <= next_si;

                    if upd_t = N_act - 1 then
                        -- Last training sample; set forecast season index
                        if si_counter = M_act_m1 then fc_si_counter <= 0;
                        else                          fc_si_counter <= si_counter + 1;
                        end if;
                        k_cnt <= 0;
                        state <= S_FORECAST;
                    else
                        -- More samples; inline prefetch (S_UPDATE_PRE skipped)
                        upd_S_pre    <= S_arr(next_si);  -- Latch next seasonal factor
                        upd_data_pre <= data_buf(t_cnt); -- Latch next sample
                        upd_si       <= next_si;
                        upd_t        <= t_cnt;
                        state        <= S_UPDATE;
                    end if;

                -- =============================================================
                -- S_FORECAST: emits one forecast per cycle for HORIZON_act steps.
                -- F(k) = L + (k+1)*T + S[fc_si_counter]
                -- =============================================================
                when S_FORECAST =>
                    k1_fp          := to_signed((k_cnt + 1) * 65536, 32); -- (k+1) in Q16.16
                    forecast_out   <= L_reg + fp_mul_data(k1_fp, T_reg)
                                           + S_arr(fc_si_counter);
                    forecast_valid <= '1';

                    if fc_si_counter = M_act_m1 then fc_si_counter <= 0;
                    else                             fc_si_counter <= fc_si_counter + 1;
                    end if;

                    if k_cnt = HORIZON_act - 1 then
                        last_forecast <= '1'; -- Final forecast step
                        state         <= S_IDLE;
                    else
                        k_cnt <= k_cnt + 1;
                    end if;

                -- =============================================================
                -- S_ERROR: pulses error_out for one cycle then returns to S_IDLE.
                -- =============================================================
                when S_ERROR =>
                    error_out <= '1';
                    state     <= S_IDLE;

                when others =>
                    state <= S_IDLE;

            end case;
        end if;
    end if;
    end process;

end architecture rtl;
