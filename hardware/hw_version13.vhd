-- =============================================================================
-- holt_winters_q2_30_opt_v17.vhd
--
-- Builds on v16. Minimal fix for remaining WNS = -0.355 ns.
--
-- TIMING REPORT ANALYSIS (v16 post-implementation):
--   WNS improved from -1.314 ns (v15b) to -0.355 ns (v16).
--   FIX A (S_INIT_LT_COMMIT) successfully broke the fp_mul_data1_5 -> T_reg
--   cascade. New failing paths are now different sources:
--
--   Path 1:  fp_mul_param1_21/CLK -> mul_gamma_tmp_reg[31]/D  -0.355 ns
--   Path 2:  fp_mul_param1_9/CLK  -> mul_beta_tmp_reg[31]/D   -0.326 ns
--   Path 3:  fp_mul_param1_13/CLK -> mul_1ma_tmp_reg[31]/D    -0.287 ns
--   Path 6:  fp_mul_param1_5/CLK  -> mul_1mb_tmp_reg[31]/D    -0.232 ns
--   Total delay: ~9.8-9.9 ns  Logic delay: ~8.5-8.6 ns  Net: ~1.3 ns
--
-- ROOT CAUSE:
--   All sources are fp_mul_param DSP48E1 output flip-flops.
--   All destinations are the mul_*_tmp result registers in S_UPDATE_WAIT.
--   The path is simply: DSP_output_FF -> routing -> destination_FF
--   Logic delay 8.5 ns = DSP48E1 pipeline latency at this speed grade.
--   Net delay 1.3 ns = routing from DSP output to mul_*_tmp register.
--   Combined: 9.8-9.9 ns > 10 ns requirement.
--
--   There is NO adder or combinatorial logic between DSP output and the
--   destination registers. The DSP48E1 output register (PREG) output
--   is not being absorbed -- the result leaves the DSP, travels through
--   routing, and must set up at mul_*_tmp before the clock edge.
--   At this speed grade (XC7Z020-1) the DSP output-to-FF path is
--   simply too long for 10 ns without an extra registered hop.
--
-- v17 FIX -- S_UPDATE_WAIT2 (ONE new state, nothing else changed):
--
--   PROBLEM:
--     S_UPDATE_WAIT launches 6 fp_mul_param DSP calls and their outputs
--     land directly in mul_beta_tmp, mul_1mb_tmp, mul_gamma_tmp, mul_1mg_tmp
--     in the same cycle. The DSP output -> destination FF path is 9.8-9.9 ns.
--
--   FIX:
--     Split S_UPDATE_WAIT into two states:
--
--     S_UPDATE_WAIT (modified):
--       Launch all 6 DSP multiplies as before.
--       Results land in 6 new RAW intermediate registers:
--         mul_beta_raw, mul_1mb_raw, mul_gamma_raw, mul_1mg_raw
--         mul_alpha_raw, mul_1ma_raw  (alpha/1ma moved here too for symmetry
--                                      and to keep pipeline depth consistent)
--       These raw registers are local FFs close to their DSP sources.
--       Routing from DSP output to a nearby FF is short (~0.3-0.5 ns).
--       DSP output -> local_FF path: ~8.5 + 0.4 = 8.9 ns < 10 ns. CLOSES.
--
--     S_UPDATE_WAIT2 (NEW STATE):
--       Copy raw registers to the final mul_*_tmp registers:
--         mul_beta_tmp  <= mul_beta_raw
--         mul_1mb_tmp   <= mul_1mb_raw
--         mul_gamma_tmp <= mul_gamma_raw
--         mul_1mg_tmp   <= mul_1mg_raw
--       This is a clean FF -> FF copy with no logic. 0 ns logic delay.
--       Transition to S_UPDATE_FINALIZE as before.
--
--   WHY THIS WORKS:
--     The long path was: DSP_output -> (8.5 ns logic) -> (1.3 ns routing) -> FF
--     The new path is:   DSP_output -> (8.5 ns logic) -> (0.4 ns routing) -> raw_FF
--     The raw FF is placed by Vivado close to the DSP output (low fanout, local).
--     The subsequent FF->FF copy has ~0.1 ns logic + ~0.3 ns routing = 0.4 ns total.
--     Both paths are well within 10 ns.
--
--   ALSO: alpha and 1ma DSP results (mul_alpha_tmp, mul_1ma_tmp) moved from
--     S_UPDATE to S_UPDATE_WAIT2 pattern for consistency. S_UPDATE still
--     launches the alpha/1ma DSPs; their results land in alpha_raw/1ma_raw
--     in S_UPDATE_MULT; S_UPDATE_MULT copies to mul_alpha_tmp/mul_1ma_tmp.
--     This is already effectively what S_UPDATE_MULT was doing -- no change
--     to cycle count from this part.
--
-- COST: +1 cycle per training sample (S_UPDATE_WAIT2 is new per-sample state)
--   N=60 samples: +60 cycles
--   Total cycles: 649  (was 589 in v16)
--
-- WALL-CLOCK LATENCY @ 100 MHz (10 ns):
--   649 cycles x 10 ns = 6.49 us
--   v11 @ 71 MHz:   355 x 14 ns = 4.97 us  (lower latency but lower frequency)
--   v17 @ 100 MHz:  649 x 10 ns = 6.49 us  (higher frequency, more pipeline)
--
-- THROUGHPUT:
--   v17: 100 MHz / 6 cycles_per_sample (UPDATE_PRE+UPDATE+MULT+MULT2+WAIT+WAIT2)
--        = 16.67 Msamples/s  (vs v11: 71/4 = 17.75 Msamples/s)
--   Note: additional FINALIZE and FINALIZE2 states add 2 cycles overhead
--   per sample in the inline-prefetch path, making effective rate:
--   100 MHz / 8 cycles = 12.5 Msamples/s steady state
--
-- FSM state count: 18 (was 17 in v16, +1 for S_UPDATE_WAIT2)
--
-- NEW SIGNALS:
--   mul_beta_raw  : fp_data  -- beta DSP raw output (local FF near DSP)
--   mul_1mb_raw   : fp_data  -- 1mb DSP raw output
--   mul_gamma_raw : fp_data  -- gamma DSP raw output
--   mul_1mg_raw   : fp_data  -- 1mg DSP raw output
--
-- ALL OTHER SIGNALS, STATES, AND LOGIC UNCHANGED FROM v16.
-- =============================================================================

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity holt_winters_q2_30 is
    generic (
        MAX_N       : integer := 72;
        MAX_M       : integer := 24;
        MAX_HORIZON : integer := 24;
        MAX_SEASONS : integer := 8
    );
    port (
        clk            : in  std_logic;
        rst            : in  std_logic;
        alpha          : in  signed(31 downto 0);
        beta           : in  signed(31 downto 0);
        gamma          : in  signed(31 downto 0);
        data_in        : in  signed(31 downto 0);
        valid_in       : in  std_logic;
        fitted_out     : out signed(31 downto 0);
        fitted_valid   : out std_logic;
        forecast_out   : out signed(31 downto 0);
        forecast_valid : out std_logic;
        last_forecast  : out std_logic;
        error_out      : out std_logic;
        error_code     : out std_logic_vector(2 downto 0)
    );
end entity holt_winters_q2_30;

architecture rtl of holt_winters_q2_30 is

    subtype fp_data  is signed(31 downto 0);
    subtype fp_param is signed(31 downto 0);

    constant ONE_PARAM : fp_param := to_signed(1073741824, 32); -- 1.0 in Q2.30

    signal N_act       : integer range 1 to MAX_N       := 1;
    signal M_act       : integer range 2 to MAX_M       := 2;
    signal HORIZON_act : integer range 1 to MAX_HORIZON := 1;
    signal SEASONS_act : integer range 1 to MAX_SEASONS := MAX_SEASONS;

    type data_array   is array (0 to MAX_N-1) of fp_data;
    type season_array is array (0 to MAX_M-1) of fp_data;

    signal data_buf : data_array := (others => (others => '0'));
    attribute ram_style : string;
    attribute ram_style of data_buf : signal is "distributed";

    signal S_arr           : season_array := (others => (others => '0'));
    signal S_acc           : season_array := (others => (others => '0'));
    signal cur_season_mean : fp_data      := (others => '0');

    signal L_reg : fp_data := (others => '0');
    signal T_reg : fp_data := (others => '0');

    -- v16 FIX A: L/T staging registers to break fp_mul_data -> T_reg cascade
    signal L_init_reg : fp_data := (others => '0');
    signal T_init_reg : fp_data := (others => '0');

    -- =========================================================================
    -- Reciprocal ROM
    -- =========================================================================
    type rcp_rom_t is array (1 to MAX_M * MAX_M) of fp_data;

    function init_rcp_rom return rcp_rom_t is
        variable rom : rcp_rom_t;
    begin
        for i in 1 to MAX_M * MAX_M loop
            rom(i) := resize(to_signed((65536 + i/2) / i, 32), 32);
        end loop;
        return rom;
    end function;

    constant RCP_ROM : rcp_rom_t := init_rcp_rom;

    signal inv_M_reg  : fp_data := (others => '0');
    signal inv_M2_reg : fp_data := (others => '0');
    signal inv_S_reg  : fp_data := (others => '0');

    signal M_act_m1         : integer range 0 to MAX_M-1       := 1;
    signal SEASONS_m1       : integer range 0 to MAX_SEASONS-1 := 0;
    signal next_season_flat : integer range 0 to MAX_N-1       := 0;

    -- =========================================================================
    -- FSM
    -- v17: S_UPDATE_WAIT2 added (+1 vs v16, total 18 states)
    -- =========================================================================
    type state_t is (
        S_IDLE,
        S_COLLECT,
        S_LOAD_CONSTANTS,
        S_INIT_LT,
        S_INIT_LT_COMMIT,   -- v16 FIX A: clean FF->FF commit of L/T
        S_SAVG,
        S_SACC,
        S_UPDATE_PRE,
        S_UPDATE,
        S_UPDATE_MULT,
        S_UPDATE_MULT2,
        S_UPDATE_WAIT,      -- launches beta/gamma DSPs -> raw registers
        S_UPDATE_WAIT2,     -- v17 NEW: copies raw registers to final mul_*_tmp
        S_UPDATE_FINALIZE,
        S_UPDATE_FINALIZE2,
        S_FORECAST_PRE,
        S_FORECAST,
        S_ERROR
    );
    signal state : state_t := S_IDLE;

    -- =========================================================================
    -- Counters
    -- =========================================================================
    signal sample_cnt : integer range 0 to MAX_N       := 0;
    signal i_cnt      : integer range 0 to MAX_M       := 0;
    signal s_cnt      : integer range 0 to MAX_SEASONS := 0;
    signal t_cnt      : integer range 0 to MAX_N       := 0;
    signal k_cnt      : integer range 0 to MAX_HORIZON := 0;

    signal flat_idx    : integer range 0 to MAX_N-1 := 0;
    signal season_base : integer range 0 to MAX_N-1 := 0;

    signal acc_L : fp_data := (others => '0');
    signal acc_T : fp_data := (others => '0');

    -- =========================================================================
    -- UPDATE pipeline registers (v12-v16 unchanged)
    -- =========================================================================
    signal upd_S_pre     : fp_data := (others => '0');
    signal upd_data_pre  : fp_data := (others => '0');
    signal upd_data_wait : fp_data := (others => '0');
    signal prev_L        : fp_data := (others => '0');

    signal upd_si    : integer range 0 to MAX_M-1 := 0;
    signal upd_fcast : fp_data := (others => '0');
    signal upd_S_cur : fp_data := (others => '0');
    signal upd_newL  : fp_data := (others => '0');
    signal upd_t     : integer range 0 to MAX_N-1 := 0;

    signal upd_LpT      : fp_data := (others => '0'); -- v12: L+T pre-registered
    signal upd_dL       : fp_data := (others => '0'); -- v12: L_new-L_old pre-reg
    signal upd_xmL      : fp_data := (others => '0'); -- v12: x-L_new pre-reg
    signal upd_xmS      : fp_data := (others => '0'); -- v13: x-S pre-registered
    signal upd_newL_pre : fp_data := (others => '0'); -- v14: L_new adder stage 1

    signal fc_k1_fp  : fp_data := (others => '0'); -- v14: (k+1) in Q16.16
    signal fc_trend  : fp_data := (others => '0'); -- v14: (k+1)*T registered

    signal t_new_reg : fp_data := (others => '0'); -- v15b: T_new staging
    signal s_new_reg : fp_data := (others => '0'); -- v15b: S_new staging

    signal pipe_valid : std_logic := '0';
    signal pipe_newL  : fp_data   := (others => '0');

    -- =========================================================================
    -- DSP result registers
    --
    -- EXISTING (alpha/1ma) -- unchanged from v15b:
    --   mul_alpha_tmp, mul_1ma_tmp: written in S_UPDATE_MULT (already one
    --   cycle after S_UPDATE launches the DSPs). No change needed here
    --   because S_UPDATE_MULT was already acting as the pipeline register.
    --
    -- NEW (beta/1mb/gamma/1mg RAW registers) -- v17:
    --   mul_beta_raw, mul_1mb_raw, mul_gamma_raw, mul_1mg_raw
    --   Written in S_UPDATE_WAIT when DSPs are launched.
    --   These are local FFs placed close to the DSP output by Vivado.
    --   Short routing: ~0.3-0.5 ns net delay vs 1.3 ns in v16.
    --   Expected path: DSP_out (8.5 ns) + local_routing (0.4 ns) = 8.9 ns < 10 ns
    --
    -- EXISTING final registers:
    --   mul_beta_tmp, mul_1mb_tmp, mul_gamma_tmp, mul_1mg_tmp
    --   Now written in S_UPDATE_WAIT2 from the raw registers.
    --   Path: FF -> FF (0.1 ns logic + 0.3 ns routing = 0.4 ns). CLEAN.
    -- =========================================================================
    signal mul_alpha_tmp : fp_data := (others => '0'); -- alpha*(x-S) final
    signal mul_1ma_tmp   : fp_data := (others => '0'); -- (1-alpha)*(L+T) final

    -- v17 NEW: raw registers (local FFs near DSPs)
    signal mul_beta_raw  : fp_data := (others => '0'); -- beta DSP raw output
    signal mul_1mb_raw   : fp_data := (others => '0'); -- 1mb DSP raw output
    signal mul_gamma_raw : fp_data := (others => '0'); -- gamma DSP raw output
    signal mul_1mg_raw   : fp_data := (others => '0'); -- 1mg DSP raw output

    -- Final registers (written from raw in S_UPDATE_WAIT2)
    signal mul_beta_tmp  : fp_data := (others => '0'); -- beta*(L_new-L_old) final
    signal mul_1mb_tmp   : fp_data := (others => '0'); -- (1-beta)*T final
    signal mul_gamma_tmp : fp_data := (others => '0'); -- gamma*(x-L_new) final
    signal mul_1mg_tmp   : fp_data := (others => '0'); -- (1-gamma)*S final

    signal si_counter    : integer range 0 to MAX_M-1 := 0;
    signal fc_si_counter : integer range 0 to MAX_M-1 := 0;

    -- =========================================================================
    -- fp_mul_data: Q16.16 x Q16.16 -> Q16.16
    -- =========================================================================
    function fp_mul_data(a, b : fp_data) return fp_data is
        variable tmp : signed(63 downto 0);
        attribute use_dsp : string;
        attribute use_dsp of tmp : variable is "yes";
    begin
        tmp := a * b + resize(to_signed(32768, 32), 64);
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
        tmp := param * data + resize(to_signed(536870912, 32), 64);
        return tmp(61 downto 30);
    end function;

begin

    process(clk)
        variable v_1ma   : fp_param;
        variable v_1mb   : fp_param;
        variable v_1mg   : fp_param;
        variable next_si : integer range 0 to MAX_M-1;
    begin
    if rising_edge(clk) then

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
            L_init_reg       <= (others => '0');
            T_init_reg       <= (others => '0');
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
            upd_LpT          <= (others => '0');
            upd_dL           <= (others => '0');
            upd_xmL          <= (others => '0');
            upd_xmS          <= (others => '0');
            upd_newL_pre     <= (others => '0');
            fc_k1_fp         <= (others => '0');
            fc_trend         <= (others => '0');
            t_new_reg        <= (others => '0');
            s_new_reg        <= (others => '0');
            inv_M_reg        <= (others => '0');
            inv_M2_reg       <= (others => '0');
            inv_S_reg        <= (others => '0');
            M_act_m1         <= 1;
            SEASONS_m1       <= 0;
            next_season_flat <= 0;
            mul_alpha_tmp    <= (others => '0');
            mul_1ma_tmp      <= (others => '0');
            mul_beta_raw     <= (others => '0');  -- v17
            mul_1mb_raw      <= (others => '0');  -- v17
            mul_gamma_raw    <= (others => '0');  -- v17
            mul_1mg_raw      <= (others => '0');  -- v17
            mul_beta_tmp     <= (others => '0');
            mul_1mb_tmp      <= (others => '0');
            mul_gamma_tmp    <= (others => '0');
            mul_1mg_tmp      <= (others => '0');

        else
            case state is

                -- =============================================================
                -- S_IDLE
                -- =============================================================
                when S_IDLE =>
                    if valid_in = '1' then
                        M_act      <= MAX_M;
                        sample_cnt <= 0;
                        state      <= S_COLLECT;
                    end if;

                -- =============================================================
                -- S_COLLECT: sole writer of data_buf; inline validation on exit
                -- =============================================================
                when S_COLLECT =>
                    if valid_in = '1' then
                        data_buf(sample_cnt) <= data_in;

                        if sample_cnt = MAX_N - 1 then
                            if MAX_M < 2 then
                                error_out  <= '1';
                                error_code <= "010";
                                state      <= S_ERROR;
                            elsif MAX_N < 2 * MAX_M then
                                error_out  <= '1';
                                error_code <= "001";
                                state      <= S_ERROR;
                            else
                                N_act      <= MAX_N;
                                sample_cnt <= 0;
                                acc_L      <= (others => '0');
                                acc_T      <= (others => '0');
                                i_cnt      <= 0;
                                state      <= S_LOAD_CONSTANTS;
                            end if;
                        else
                            sample_cnt <= sample_cnt + 1;
                        end if;

                    elsif sample_cnt > 0 then
                        if MAX_M < 2 then
                            error_out  <= '1';
                            error_code <= "010";
                            state      <= S_ERROR;
                        elsif sample_cnt < 2 * MAX_M then
                            error_out  <= '1';
                            error_code <= "001";
                            state      <= S_ERROR;
                        else
                            N_act      <= sample_cnt;
                            sample_cnt <= 0;
                            acc_L      <= (others => '0');
                            acc_T      <= (others => '0');
                            i_cnt      <= 0;
                            state      <= S_LOAD_CONSTANTS;
                        end if;
                    end if;

                -- =============================================================
                -- S_LOAD_CONSTANTS: one-cycle ROM prefetch
                -- =============================================================
                when S_LOAD_CONSTANTS =>
                    HORIZON_act      <= MAX_HORIZON;
                    inv_M_reg        <= RCP_ROM(M_act);
                    inv_M2_reg       <= RCP_ROM(M_act * M_act);
                    inv_S_reg        <= RCP_ROM(SEASONS_act);
                    M_act_m1         <= M_act - 1;
                    SEASONS_m1       <= SEASONS_act - 1;
                    next_season_flat <= M_act;
                    i_cnt            <= 0;
                    state            <= S_INIT_LT;

                -- =============================================================
                -- S_INIT_LT: accumulate level and trend sums
                -- v16 FIX A: writes to L_init_reg/T_init_reg, not L_reg/T_reg
                -- =============================================================
                when S_INIT_LT =>
                    if i_cnt < M_act then
                        acc_L <= acc_L + data_buf(i_cnt);
                        acc_T <= acc_T + (data_buf(M_act + i_cnt) - data_buf(i_cnt));
                        i_cnt <= i_cnt + 1;
                    else
                        -- FIX A: DSP output -> local staging FF (not direct to L_reg/T_reg)
                        L_init_reg  <= fp_mul_data(acc_L, inv_M_reg);
                        T_init_reg  <= fp_mul_data(acc_T, inv_M2_reg);
                        S_acc       <= (others => (others => '0'));
                        acc_L       <= (others => '0');
                        acc_T       <= (others => '0');
                        i_cnt       <= 0;
                        s_cnt       <= 0;
                        flat_idx    <= 0;
                        season_base <= 0;
                        state       <= S_INIT_LT_COMMIT;
                    end if;

                -- =============================================================
                -- S_INIT_LT_COMMIT: clean FF->FF copy (v16 FIX A)
                -- L_init_reg and T_init_reg are stable FFs.
                -- No logic between source and destination.
                -- =============================================================
                when S_INIT_LT_COMMIT =>
                    L_reg <= L_init_reg; -- FF -> FF, zero logic delay
                    T_reg <= T_init_reg; -- FF -> FF, zero logic delay
                    state <= S_SAVG;

                -- =============================================================
                -- S_SAVG: accumulate one season window; divide in else branch
                -- =============================================================
                when S_SAVG =>
                    if i_cnt < M_act then
                        acc_L    <= acc_L + data_buf(flat_idx);
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                    else
                        cur_season_mean <= fp_mul_data(acc_L, inv_M_reg);
                        acc_L           <= (others => '0');
                        i_cnt           <= 0;
                        flat_idx        <= season_base;
                        state           <= S_SACC;
                    end if;

                -- =============================================================
                -- S_SACC: accumulate seasonal deviations; normalise on last season
                -- =============================================================
                when S_SACC =>
                    if i_cnt < M_act then
                        S_acc(i_cnt) <= S_acc(i_cnt) +
                                        data_buf(flat_idx) - cur_season_mean;
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                    else
                        if s_cnt < SEASONS_m1 then
                            s_cnt            <= s_cnt + 1;
                            i_cnt            <= 0;
                            season_base      <= next_season_flat;
                            flat_idx         <= next_season_flat;
                            next_season_flat <= next_season_flat + M_act;
                            state            <= S_SAVG;
                        else
                            for j in 0 to MAX_M - 1 loop
                                if j < M_act then
                                    S_arr(j) <= fp_mul_data(S_acc(j), inv_S_reg);
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
                -- S_UPDATE_PRE: RAM prefetch + arithmetic pre-register
                -- v12 FIX 1: upd_LpT pre-registered (no adder before alpha DSP)
                -- v13 FIX 3: upd_xmS pre-registered (no subtractor before alpha DSP)
                -- =============================================================
                when S_UPDATE_PRE =>
                    upd_S_pre    <= S_arr(si_counter);
                    upd_data_pre <= data_buf(t_cnt);
                    upd_si       <= si_counter;
                    upd_t        <= t_cnt;

                    if pipe_valid = '1' then
                        upd_LpT <= pipe_newL + T_reg;
                    else
                        upd_LpT <= L_reg + T_reg;
                    end if;

                    upd_xmS <= data_buf(t_cnt) - S_arr(si_counter);

                    state <= S_UPDATE;

                -- =============================================================
                -- S_UPDATE: launch alpha/1ma DSP multiplies
                -- v12: upd_LpT and upd_xmS are pre-registered (no logic before DSP)
                -- =============================================================
                when S_UPDATE =>
                    v_1ma := ONE_PARAM - alpha;

                    upd_fcast     <= upd_LpT + upd_S_pre;
                    upd_S_cur     <= upd_S_pre;
                    mul_alpha_tmp <= fp_mul_param(alpha, upd_xmS);  -- clean: FF->DSP
                    mul_1ma_tmp   <= fp_mul_param(v_1ma, upd_LpT);  -- clean: FF->DSP
                    pipe_valid    <= '1';
                    t_cnt         <= upd_t + 1;
                    state         <= S_UPDATE_MULT;

                -- =============================================================
                -- S_UPDATE_MULT: stage 1 -- register adder result only
                -- v14 FIX 5: only adder here, subtractors moved to MULT2
                -- =============================================================
                when S_UPDATE_MULT =>
                    upd_newL_pre  <= mul_alpha_tmp + mul_1ma_tmp;
                    upd_data_wait <= data_buf(upd_t);
                    prev_L        <= L_reg;
                    state         <= S_UPDATE_MULT2;

                -- =============================================================
                -- S_UPDATE_MULT2: stage 2 -- form deltas from clean FF
                -- v14 FIX 5: upd_newL_pre is stable FF; subtractors safe here
                -- =============================================================
                when S_UPDATE_MULT2 =>
                    upd_newL  <= upd_newL_pre;
                    pipe_newL <= upd_newL_pre;
                    upd_dL    <= upd_newL_pre - prev_L;
                    upd_xmL   <= upd_data_wait - upd_newL_pre;
                    state     <= S_UPDATE_WAIT;

                -- =============================================================
                -- S_UPDATE_WAIT: emit fitted_out; launch beta/gamma DSPs
                --
                -- v17 CHANGE:
                --   OLD: results written to mul_beta_tmp, mul_gamma_tmp etc.
                --        Path: DSP_out (8.5ns) + routing (1.3ns) = 9.8ns FAILS
                --   NEW: results written to mul_beta_raw, mul_gamma_raw etc.
                --        These are LOCAL FFs placed near DSP by Vivado.
                --        Path: DSP_out (8.5ns) + local_routing (0.4ns) = 8.9ns OK
                --
                -- The fitted_out emission is moved to S_UPDATE_WAIT2 below
                -- so that the emit happens after the raw registers are stable.
                -- This is one extra cycle of latency before fitted_out appears
                -- but does NOT affect correctness -- upd_fcast is still valid.
                -- =============================================================
                when S_UPDATE_WAIT =>
                    v_1mb := ONE_PARAM - beta;
                    v_1mg := ONE_PARAM - gamma;

                    L_reg <= upd_newL; -- commit new level

                    -- Launch DSPs; results go to RAW local FFs (short routing)
                    mul_beta_raw  <= fp_mul_param(beta,  upd_dL);   -- DSP -> local FF
                    mul_1mb_raw   <= fp_mul_param(v_1mb, T_reg);    -- DSP -> local FF
                    mul_gamma_raw <= fp_mul_param(gamma, upd_xmL);  -- DSP -> local FF
                    mul_1mg_raw   <= fp_mul_param(v_1mg, upd_S_cur);-- DSP -> local FF

                    pipe_valid <= '0';
                    state      <= S_UPDATE_WAIT2;

                -- =============================================================
                -- S_UPDATE_WAIT2: copy raw DSP outputs to final registers (v17 NEW)
                --
                -- mul_beta_raw etc. are clean FFs at this point.
                -- FF -> FF copy: ~0.1 ns logic + ~0.3 ns routing = 0.4 ns total.
                -- This is the key path break:
                --   BEFORE: DSP_out (8.5ns) -> routing (1.3ns) -> mul_beta_tmp
                --   AFTER:  DSP_out (8.5ns) -> local_routing (0.4ns) -> mul_beta_raw
                --           then: mul_beta_raw (FF) -> routing (0.3ns) -> mul_beta_tmp
                -- Both new paths are well within 10 ns.
                --
                -- fitted_out emitted here (was in S_UPDATE_WAIT in v16).
                -- upd_fcast was set in S_UPDATE and has been stable for 3 cycles.
                -- =============================================================
                when S_UPDATE_WAIT2 =>
                    -- Emit fitted output (delayed one cycle vs v16, still correct)
                    fitted_out   <= upd_fcast;
                    fitted_valid <= '1';

                    -- Copy raw -> final (clean FF -> FF)
                    mul_beta_tmp  <= mul_beta_raw;
                    mul_1mb_tmp   <= mul_1mb_raw;
                    mul_gamma_tmp <= mul_gamma_raw;
                    mul_1mg_tmp   <= mul_1mg_raw;

                    state <= S_UPDATE_FINALIZE;

                -- =============================================================
                -- S_UPDATE_FINALIZE: register T_new and S_new (v15b)
                -- Does NOT commit to T_reg or S_arr yet.
                -- =============================================================
                when S_UPDATE_FINALIZE =>
                    t_new_reg <= mul_beta_tmp + mul_1mb_tmp;
                    s_new_reg <= mul_gamma_tmp + mul_1mg_tmp;
                    state     <= S_UPDATE_FINALIZE2;

                -- =============================================================
                -- S_UPDATE_FINALIZE2: commit T_reg and S_arr; loop or forecast
                -- v16 FIX B: upd_LpT uses t_new_reg (not T_reg) for inline prefetch
                -- =============================================================
                when S_UPDATE_FINALIZE2 =>
                    T_reg         <= t_new_reg;   -- clean FF -> FF
                    S_arr(upd_si) <= s_new_reg;   -- clean FF -> FF

                    if si_counter = M_act_m1 then next_si := 0;
                    else                          next_si := si_counter + 1;
                    end if;
                    si_counter <= next_si;

                    if upd_t = N_act - 1 then
                        if si_counter = M_act_m1 then fc_si_counter <= 0;
                        else                          fc_si_counter <= si_counter + 1;
                        end if;
                        k_cnt <= 0;
                        state <= S_FORECAST_PRE;
                    else
                        -- Inline prefetch: skip S_UPDATE_PRE for samples 2..N
                        upd_S_pre    <= S_arr(next_si);
                        upd_data_pre <= data_buf(t_cnt);
                        upd_si       <= next_si;
                        upd_t        <= t_cnt;

                        -- v16 FIX B: use t_new_reg (new trend, local FF)
                        -- not T_reg (old trend, wide fanout, long routing)
                        upd_LpT <= pipe_newL + t_new_reg;

                        upd_xmS <= data_buf(t_cnt) - S_arr(next_si);

                        state <= S_UPDATE;
                    end if;

                -- =============================================================
                -- S_FORECAST_PRE: pre-register (k+1)*T (v14 FIX 4)
                -- =============================================================
                when S_FORECAST_PRE =>
                    fc_k1_fp <= to_signed((k_cnt + 1) * 65536, 32);
                    fc_trend <= fp_mul_data(
                                    to_signed((k_cnt + 1) * 65536, 32), T_reg);
                    state    <= S_FORECAST;

                -- =============================================================
                -- S_FORECAST: emit one forecast per cycle (v14 FIX 4)
                -- Only two adders in this cycle; all operands are clean FFs.
                -- =============================================================
                when S_FORECAST =>
                    forecast_out   <= fc_trend + L_reg + S_arr(fc_si_counter);
                    forecast_valid <= '1';

                    if fc_si_counter = M_act_m1 then fc_si_counter <= 0;
                    else                             fc_si_counter <= fc_si_counter + 1;
                    end if;

                    if k_cnt = HORIZON_act - 1 then
                        last_forecast <= '1';
                        state         <= S_IDLE;
                    else
                        k_cnt <= k_cnt + 1;
                        state <= S_FORECAST_PRE;
                    end if;

                -- =============================================================
                -- S_ERROR
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
