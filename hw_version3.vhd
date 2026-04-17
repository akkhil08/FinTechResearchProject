-- =============================================================================
-- holt_winters_q2_30_fixed.vhd
--
-- FIXES APPLIED:
--
-- FIX 1 -- Serialize data_buf reads (reduces RAM64M from 132 to 11)
--   S_INIT_LT + S_INIT_LT_ODD : max 2 unique reads/cycle
--   S_SAVG    + S_SAVG_ODD    : max 1 unique read/cycle
--   S_SACC    + S_SACC_ODD    : max 1 unique read/cycle
--   Global worst case: 2 unique reads -> 1 copy sufficient
--
-- FIX 2 -- DSP48 on architecture-level signals
--   Vivado ignores use_dsp on variables inside functions.
--   Six architecture-level signals with use_dsp="yes" correctly
--   infer DSP48E1. Each saves ~7ns vs LUT multiply.
--
-- FIX 3 -- Pipeline stages for S_UPDATE
--   S_UPDATE      : latch side signals, read data_buf (1 read)
--   S_UPDATE_MUL  : LAUNCH DSP multiplications only
--   S_UPDATE_MUL2 : CAPTURE DSP results (fixes one-cycle-behind bug)
--   S_UPDATE_WAIT : T_reg and S_arr updates
--
-- FIX 3 KEY BUG FIXED -- DSP one-cycle-behind error
--   PROBLEM: S_UPDATE_MUL both launched AND read dsp_alpha_x_data
--            in the same cycle. On first call dsp signals are zero
--            from reset, so upd_newL = 0 for sample 0.
--            Result: L_reg was always computed from wrong data,
--            causing forecast values of ~7360 instead of ~7062.
--   SOLUTION: S_UPDATE_MUL  -> launches multiply only
--             S_UPDATE_MUL2 -> captures result (fully registered)
--             Two states means DSP output is always valid when read.
--
-- FIX 4 -- Eliminate second data_buf read in S_UPDATE_WAIT
--   gamma path reconstructed as upd1_data + upd1_S_cur - upd_newL
--   avoids a second RAM read entirely.
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
        clk             : in  std_logic;
        rst             : in  std_logic;
        m_in            : in  std_logic_vector(4 downto 0);
        forecast_start  : in  std_logic;
        horizon_in      : in  std_logic_vector(4 downto 0);
        start           : in  std_logic;
        alpha           : in  signed(31 downto 0);
        beta            : in  signed(31 downto 0);
        gamma           : in  signed(31 downto 0);
        data_in         : in  signed(31 downto 0);
        valid_in        : in  std_logic;
        fitted_out      : out signed(31 downto 0);
        fitted_valid    : out std_logic;
        forecast_out    : out signed(31 downto 0);
        forecast_valid  : out std_logic;
        last_forecast   : out std_logic;
        error_out       : out std_logic;
        error_code      : out std_logic_vector(2 downto 0)
    );
end entity holt_winters_q2_30;

architecture rtl of holt_winters_q2_30 is

    subtype fp_data  is signed(31 downto 0);
    subtype fp_param is signed(31 downto 0);

    constant ONE_PARAM : fp_param := to_signed(1073741824, 32);

    signal N_act       : integer range 1 to MAX_N       := 1;
    signal M_act       : integer range 2 to MAX_M       := 2;
    signal HORIZON_act : integer range 1 to MAX_HORIZON := 1;
    signal SEASONS_act : integer range 1 to MAX_SEASONS := 1;

    type data_array   is array (0 to MAX_N-1) of fp_data;
    type season_array is array (0 to MAX_M-1) of fp_data;

    signal data_buf        : data_array   := (others => (others => '0'));
    signal S_arr           : season_array := (others => (others => '0'));
    signal S_acc           : season_array := (others => (others => '0'));
    signal S_acc2          : season_array := (others => (others => '0'));
    signal cur_season_mean : fp_data      := (others => '0');
    signal L_reg           : fp_data      := (others => '0');
    signal T_reg           : fp_data      := (others => '0');

    -- =========================================================================
    -- FSM
    -- FIX 1: Added ODD states for serialized reads
    -- FIX 3: Added S_UPDATE_MUL2 to capture DSP results one cycle after launch
    -- =========================================================================
    type state_t is (
        S_IDLE,
        S_COLLECT,
        S_VALIDATE,
        S_WAIT_FC,
        S_INIT_LT,      -- reads even: data_buf(i_cnt), data_buf(M+i_cnt)
        S_INIT_LT_ODD,  -- reads odd:  data_buf(i_cnt+1), data_buf(M+i_cnt+1)
        S_INIT_LT_DIV,  -- computes L0, T0
        S_SAVG,         -- reads even: data_buf(idx)
        S_SAVG_ODD,     -- reads odd:  data_buf(idx+1)
        S_SAVG_DIV,     -- computes season mean
        S_SACC,         -- reads even: data_buf(idx)
        S_SACC_ODD,     -- reads odd:  data_buf(idx+1)
        S_UPDATE,       -- reads data_buf(t_cnt), latches side signals
        S_UPDATE_MUL,   -- LAUNCHES DSP multiplications (no read of results)
        S_UPDATE_MUL2,  -- CAPTURES DSP results (fully registered, correct)
        S_UPDATE_WAIT,  -- T_reg and S_arr updates using captured DSP results
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

    -- =========================================================================
    -- Accumulators
    -- =========================================================================
    signal acc_L  : fp_data := (others => '0');
    signal acc_T  : fp_data := (others => '0');
    signal acc_L2 : fp_data := (others => '0');
    signal acc_T2 : fp_data := (others => '0');

    -- =========================================================================
    -- UPDATE pipeline registers
    --
    -- Stage 1 (S_UPDATE):
    --   upd1_data  = data_buf(t_cnt) - S_arr(si)
    --   upd1_LpT   = L+T or pipe_newL+T
    --   upd1_S_cur = S_arr(si)
    --   upd1_fcast = LpT + S_arr(si)
    --   upd1_si    = si_counter value
    --   upd1_t     = t_cnt value
    --
    -- Stage 2 (S_UPDATE_MUL): launches DSP, passes upd2_* through
    --
    -- Stage 3 (S_UPDATE_MUL2): captures DSP results into upd_newL
    -- =========================================================================
    signal upd1_data  : fp_data := (others => '0');
    signal upd1_LpT   : fp_data := (others => '0');
    signal upd1_S_cur : fp_data := (others => '0');
    signal upd1_fcast : fp_data := (others => '0');
    signal upd1_si    : integer range 0 to MAX_M-1 := 0;
    signal upd1_t     : integer range 0 to MAX_N-1 := 0;

    -- passed through S_UPDATE_MUL unchanged
    signal upd2_data  : fp_data := (others => '0');
    signal upd2_LpT   : fp_data := (others => '0');
    signal upd2_S_cur : fp_data := (others => '0');
    signal upd2_fcast : fp_data := (others => '0');
    signal upd2_si    : integer range 0 to MAX_M-1 := 0;
    signal upd2_t     : integer range 0 to MAX_N-1 := 0;

    -- captured in S_UPDATE_MUL2 from DSP outputs
    signal upd_si    : integer range 0 to MAX_M-1 := 0;
    signal upd_fcast : fp_data := (others => '0');
    signal upd_S_cur : fp_data := (others => '0');
    signal upd_newL  : fp_data := (others => '0');
    signal upd_t     : integer range 0 to MAX_N-1 := 0;
    signal upd_data  : fp_data := (others => '0');

    signal pipe_valid : std_logic := '0';
    signal pipe_newL  : fp_data   := (others => '0');

    -- =========================================================================
    -- Season index counters
    -- =========================================================================
    signal si_counter    : integer range 0 to MAX_M-1 := 0;
    signal fc_si_counter : integer range 0 to MAX_M-1 := 0;

    -- =========================================================================
    -- FIX 2: DSP48 signals at architecture level
    -- Six signals covering all multiply operations on the critical path.
    -- use_dsp on architecture signals is correctly honoured by Vivado.
    -- use_dsp on variables inside functions is silently ignored.
    --
    -- S_UPDATE_MUL  launches: dsp_alpha_x_data, dsp_1ma_x_LpT
    -- S_UPDATE_MUL2 launches: dsp_beta_x_dL, dsp_1mb_x_T,
    --                         dsp_gamma_x_dS, dsp_1mg_x_Scur
    -- S_UPDATE_WAIT reads results of S_UPDATE_MUL2 launches
    -- =========================================================================
    attribute use_dsp : string;

    signal dsp_alpha_x_data : signed(63 downto 0) := (others => '0');
    signal dsp_1ma_x_LpT    : signed(63 downto 0) := (others => '0');
    signal dsp_beta_x_dL    : signed(63 downto 0) := (others => '0');
    signal dsp_1mb_x_T      : signed(63 downto 0) := (others => '0');
    signal dsp_gamma_x_dS   : signed(63 downto 0) := (others => '0');
    signal dsp_1mg_x_Scur   : signed(63 downto 0) := (others => '0');

    attribute use_dsp of dsp_alpha_x_data : signal is "yes";
    attribute use_dsp of dsp_1ma_x_LpT   : signal is "yes";
    attribute use_dsp of dsp_beta_x_dL   : signal is "yes";
    attribute use_dsp of dsp_1mb_x_T     : signal is "yes";
    attribute use_dsp of dsp_gamma_x_dS  : signal is "yes";
    attribute use_dsp of dsp_1mg_x_Scur  : signal is "yes";

    -- =========================================================================
    -- Helper functions (init states only -- not on critical path)
    -- =========================================================================
    function fp_mul_data(a, b : fp_data) return fp_data is
        variable tmp : signed(63 downto 0);
    begin
        tmp := a * b + resize(to_signed(32768, 32), 64);
        return tmp(47 downto 16);
    end function;

    function fp_mul_param(param : fp_param; data : fp_data) return fp_data is
        variable tmp : signed(63 downto 0);
    begin
        tmp := param * data + resize(to_signed(536870912, 32), 64);
        return tmp(61 downto 30);
    end function;

    function rcp(x : integer) return fp_data is
    begin
        return resize(to_signed((65536 + x/2) / x, 32), 32);
    end function;

    function divides(a, b : integer) return boolean is
    begin
        return (b mod a) = 0;
    end function;

begin

    process(clk)
        variable v_1ma       : fp_param;
        variable v_1mb       : fp_param;
        variable v_1mg       : fp_param;
        variable v_si        : integer range 0 to MAX_M-1;
        variable inv_M       : fp_data;
        variable inv_M2      : fp_data;
        variable inv_SEASONS : fp_data;
        variable idx         : integer;
        variable k1_fp       : fp_data;
        variable v_mean      : fp_data;
    begin
    if rising_edge(clk) then

        forecast_valid <= '0';
        fitted_valid   <= '0';
        last_forecast  <= '0';
        error_out      <= '0';
        error_code     <= (others => '0');

        if rst = '1' then
            state           <= S_IDLE;
            sample_cnt      <= 0; i_cnt <= 0; s_cnt <= 0;
            t_cnt           <= 0; k_cnt <= 0;
            acc_L           <= (others => '0');
            acc_L2          <= (others => '0');
            acc_T           <= (others => '0');
            acc_T2          <= (others => '0');
            L_reg           <= (others => '0');
            T_reg           <= (others => '0');
            pipe_valid      <= '0';
            pipe_newL       <= (others => '0');
            S_acc           <= (others => (others => '0'));
            S_acc2          <= (others => (others => '0'));
            S_arr           <= (others => (others => '0'));
            cur_season_mean <= (others => '0');
            si_counter      <= 0;
            fc_si_counter   <= 0;
            upd1_data       <= (others => '0');
            upd1_LpT        <= (others => '0');
            upd1_S_cur      <= (others => '0');
            upd1_fcast      <= (others => '0');
            upd1_si         <= 0;
            upd1_t          <= 0;
            upd2_data       <= (others => '0');
            upd2_LpT        <= (others => '0');
            upd2_S_cur      <= (others => '0');
            upd2_fcast      <= (others => '0');
            upd2_si         <= 0;
            upd2_t          <= 0;
            upd_si          <= 0;
            upd_fcast       <= (others => '0');
            upd_S_cur       <= (others => '0');
            upd_newL        <= (others => '0');
            upd_t           <= 0;
            upd_data        <= (others => '0');
            dsp_alpha_x_data <= (others => '0');
            dsp_1ma_x_LpT   <= (others => '0');
            dsp_beta_x_dL   <= (others => '0');
            dsp_1mb_x_T     <= (others => '0');
            dsp_gamma_x_dS  <= (others => '0');
            dsp_1mg_x_Scur  <= (others => '0');

        else
            case state is

                -- =============================================================
                when S_IDLE =>
                    if start = '1' then
                        M_act      <= to_integer(unsigned(m_in));
                        sample_cnt <= 0;
                        state      <= S_COLLECT;
                    end if;

                -- =============================================================
                when S_COLLECT =>
                    if valid_in = '1' then
                        data_buf(sample_cnt) <= data_in;
                        if sample_cnt = MAX_N - 1 then
                            N_act <= MAX_N; sample_cnt <= 0;
                            acc_L <= (others => '0'); acc_L2 <= (others => '0');
                            acc_T <= (others => '0'); acc_T2 <= (others => '0');
                            i_cnt <= 0; state <= S_VALIDATE;
                        else
                            sample_cnt <= sample_cnt + 1;
                        end if;
                    elsif sample_cnt > 0 then
                        N_act <= sample_cnt; sample_cnt <= 0;
                        acc_L <= (others => '0'); acc_L2 <= (others => '0');
                        acc_T <= (others => '0'); acc_T2 <= (others => '0');
                        i_cnt <= 0; state <= S_VALIDATE;
                    end if;

                -- =============================================================
                when S_VALIDATE =>
                    if    M_act < 2                     then
                        error_out <= '1'; error_code <= "010"; state <= S_ERROR;
                    elsif N_act < 2 * M_act             then
                        error_out <= '1'; error_code <= "001"; state <= S_ERROR;
                    elsif not divides(M_act, N_act)     then
                        error_out <= '1'; error_code <= "011"; state <= S_ERROR;
                    elsif (N_act / M_act) > MAX_SEASONS then
                        error_out <= '1'; error_code <= "101"; state <= S_ERROR;
                    else
                        SEASONS_act <= N_act / M_act;
                        state <= S_WAIT_FC;
                    end if;

                -- =============================================================
                when S_WAIT_FC =>
                    if forecast_start = '1' then
                        HORIZON_act <= to_integer(unsigned(horizon_in));
                        state       <= S_INIT_LT;
                    end if;

                -- =============================================================
                -- S_INIT_LT: reads EVEN index pair only
                -- Reads: data_buf(i_cnt) and data_buf(M+i_cnt) = 2 unique reads
                -- data_buf(i_cnt) reused for subtraction = not a new port
                -- Goes to S_INIT_LT_ODD for odd index pair
                -- Goes to S_INIT_LT_DIV when i_cnt >= M_act
                -- =============================================================
                when S_INIT_LT =>
                    if i_cnt < M_act then
                        acc_L <= acc_L + data_buf(i_cnt);
                        acc_T <= acc_T + (data_buf(M_act + i_cnt)
                                        - data_buf(i_cnt));
                        state <= S_INIT_LT_ODD;
                    else
                        -- All pairs done -- combine and go to divide
                        acc_L  <= acc_L + acc_L2;
                        acc_T  <= acc_T + acc_T2;
                        acc_L2 <= (others => '0');
                        acc_T2 <= (others => '0');
                        i_cnt  <= 0;
                        state  <= S_INIT_LT_DIV;
                    end if;

                -- =============================================================
                -- S_INIT_LT_ODD: reads ODD index pair only
                -- Reads: data_buf(i_cnt+1) and data_buf(M+i_cnt+1) = 2 unique
                -- Skips if i_cnt+1 >= M_act (handles odd M)
                -- Advances i_cnt by 2 then returns to S_INIT_LT
                -- =============================================================
                when S_INIT_LT_ODD =>
                    if i_cnt + 1 < M_act then
                        acc_L2 <= acc_L2 + data_buf(i_cnt + 1);
                        acc_T2 <= acc_T2 + (data_buf(M_act + i_cnt + 1)
                                          - data_buf(i_cnt + 1));
                    end if;
                    i_cnt <= i_cnt + 2;
                    state <= S_INIT_LT;

                -- =============================================================
                -- S_INIT_LT_DIV: compute L0 and T0 -- no data_buf reads
                -- =============================================================
                when S_INIT_LT_DIV =>
                    inv_M  := rcp(M_act);
                    inv_M2 := rcp(M_act * M_act);
                    L_reg  <= fp_mul_data(acc_L, inv_M);
                    T_reg  <= fp_mul_data(acc_T, inv_M2);
                    S_acc  <= (others => (others => '0'));
                    S_acc2 <= (others => (others => '0'));
                    acc_L  <= (others => '0'); acc_L2 <= (others => '0');
                    acc_T  <= (others => '0'); acc_T2 <= (others => '0');
                    i_cnt  <= 0; s_cnt <= 0;
                    state  <= S_SAVG;

                -- =============================================================
                -- S_SAVG: reads EVEN index only = 1 unique read
                -- Goes to S_SAVG_ODD for odd index
                -- Goes to S_SAVG_DIV when all elements done
                -- =============================================================
                when S_SAVG =>
                    idx := s_cnt * M_act + i_cnt;
                    if i_cnt < M_act then
                        acc_L <= acc_L + data_buf(idx);
                        state <= S_SAVG_ODD;
                    else
                        acc_L  <= acc_L + acc_L2;
                        acc_L2 <= (others => '0');
                        i_cnt  <= 0;
                        state  <= S_SAVG_DIV;
                    end if;

                -- =============================================================
                -- S_SAVG_ODD: reads ODD index only = 1 unique read
                -- Skips if i_cnt+1 >= M_act (handles odd M)
                -- Advances i_cnt by 2 then returns to S_SAVG
                -- =============================================================
                when S_SAVG_ODD =>
                    idx := s_cnt * M_act + i_cnt + 1;
                    if i_cnt + 1 < M_act then
                        acc_L2 <= acc_L2 + data_buf(idx);
                    end if;
                    i_cnt <= i_cnt + 2;
                    state <= S_SAVG;

                -- =============================================================
                -- S_SAVG_DIV: compute season mean -- no data_buf reads
                -- =============================================================
                when S_SAVG_DIV =>
                    inv_M           := rcp(M_act);
                    v_mean          := fp_mul_data(acc_L, inv_M);
                    cur_season_mean <= v_mean;
                    acc_L           <= (others => '0');
                    acc_L2          <= (others => '0');
                    i_cnt           <= 0;
                    state           <= S_SACC;

                -- =============================================================
                -- S_SACC: reads EVEN index only = 1 unique read
                -- Goes to S_SACC_ODD for odd index
                -- Finalises S_arr on last season then goes to S_UPDATE
                -- =============================================================
                when S_SACC =>
                    idx := s_cnt * M_act + i_cnt;
                    if i_cnt < M_act then
                        S_acc(i_cnt) <= S_acc(i_cnt)
                                      + data_buf(idx)
                                      - cur_season_mean;
                        state <= S_SACC_ODD;
                    else
                        -- All elements for this season done
                        if s_cnt < SEASONS_act - 1 then
                            s_cnt <= s_cnt + 1;
                            i_cnt <= 0;
                            state <= S_SAVG;
                        else
                            -- Final season -- compute S_arr
                            inv_SEASONS := rcp(SEASONS_act);
                            for j in 0 to MAX_M - 1 loop
                                if j < M_act then
                                    S_arr(j) <= fp_mul_data(
                                        S_acc(j) + S_acc2(j), inv_SEASONS);
                                end if;
                            end loop;
                            si_counter    <= 0;
                            fc_si_counter <= 0;
                            i_cnt         <= 0;
                            s_cnt         <= 0;
                            t_cnt         <= 0;
                            pipe_valid    <= '0';
                            state         <= S_UPDATE;
                        end if;
                    end if;

                -- =============================================================
                -- S_SACC_ODD: reads ODD index only = 1 unique read
                -- Skips if i_cnt+1 >= M_act (handles odd M)
                -- Advances i_cnt by 2 then returns to S_SACC
                -- =============================================================
                when S_SACC_ODD =>
                    idx := s_cnt * M_act + i_cnt + 1;
                    if i_cnt + 1 < M_act then
                        S_acc2(i_cnt + 1) <= S_acc2(i_cnt + 1)
                                           + data_buf(idx)
                                           - cur_season_mean;
                    end if;
                    i_cnt <= i_cnt + 2;
                    state <= S_SACC;

                -- =============================================================
                -- S_UPDATE: PIPELINE STAGE 1
                --
                -- Reads data_buf(t_cnt) -- 1 unique read only.
                -- Latches all side-channel values into upd1_* registers.
                -- No multiply here.
                --
                -- Path: RAM read(2ns) + subtract(1ns) + add(1ns) = 4ns  OK
                -- =============================================================
                when S_UPDATE =>
                    v_si := si_counter;

                    -- Single read from data_buf
                    upd1_data  <= data_buf(t_cnt) - S_arr(v_si);

                    -- Latch side signals -- no RAM read needed
                    upd1_S_cur <= S_arr(v_si);
                    upd1_si    <= v_si;
                    upd1_t     <= t_cnt;

                    -- Simple adds only
                    if pipe_valid = '1' then
                        upd1_LpT   <= pipe_newL + T_reg;
                        upd1_fcast <= pipe_newL + T_reg + S_arr(v_si);
                    else
                        upd1_LpT   <= L_reg + T_reg;
                        upd1_fcast <= L_reg + T_reg + S_arr(v_si);
                    end if;

                    t_cnt <= t_cnt + 1;

                    if si_counter = M_act - 1 then
                        si_counter <= 0;
                    else
                        si_counter <= si_counter + 1;
                    end if;

                    state <= S_UPDATE_MUL;

                -- =============================================================
                -- S_UPDATE_MUL: PIPELINE STAGE 2
                --
                -- LAUNCHES DSP multiplications only.
                -- Does NOT read the DSP results -- they are not ready yet.
                -- Results are captured in S_UPDATE_MUL2 next cycle.
                --
                -- This separation is the key fix for the one-cycle-behind bug:
                --   BEFORE: launched and read in same cycle -> read stale zero
                --   AFTER:  launch here, read in MUL2      -> always correct
                --
                -- Path: register inputs + launch DSP = 2ns  OK
                -- =============================================================
                when S_UPDATE_MUL =>
                    v_1ma := ONE_PARAM - alpha;

                    -- LAUNCH multiplications -- results ready NEXT cycle
                    dsp_alpha_x_data <= alpha * upd1_data
                                      + resize(to_signed(536870912, 32), 64);
                    dsp_1ma_x_LpT   <= v_1ma * upd1_LpT
                                      + resize(to_signed(536870912, 32), 64);

                    -- Pass stage 1 values through to stage 3
                    upd2_data  <= upd1_data;
                    upd2_LpT   <= upd1_LpT;
                    upd2_S_cur <= upd1_S_cur;
                    upd2_fcast <= upd1_fcast;
                    upd2_si    <= upd1_si;
                    upd2_t     <= upd1_t;

                    state <= S_UPDATE_MUL2;

                -- =============================================================
                -- S_UPDATE_MUL2: PIPELINE STAGE 3 -- NEW STATE
                --
                -- CAPTURES DSP results from S_UPDATE_MUL.
                -- dsp_alpha_x_data and dsp_1ma_x_LpT are now fully registered
                -- and hold the correct values for this sample.
                --
                -- Also launches DSP multiplications for T and S updates
                -- so S_UPDATE_WAIT only needs to read registered results.
                --
                -- Path: read registered DSP + add + launch next DSP = 5ns  OK
                -- =============================================================
                when S_UPDATE_MUL2 =>
                    v_1mb := ONE_PARAM - beta;
                    v_1mg := ONE_PARAM - gamma;

                    -- Capture newL from now-valid DSP results
                    upd_newL  <= dsp_alpha_x_data(61 downto 30)
                               + dsp_1ma_x_LpT(61 downto 30);

                    -- Update pipe_newL for next iteration
                    pipe_newL  <= dsp_alpha_x_data(61 downto 30)
                                + dsp_1ma_x_LpT(61 downto 30);
                    pipe_valid <= '1';

                    -- Launch T and S DSP multiplications
                    -- upd_newL not yet registered so use DSP slices directly
                    -- for the dL term: newL - L_reg
                    dsp_beta_x_dL  <= beta  * (dsp_alpha_x_data(61 downto 30)
                                             + dsp_1ma_x_LpT(61 downto 30)
                                             - L_reg)
                                    + resize(to_signed(536870912, 32), 64);
                    dsp_1mb_x_T    <= v_1mb * T_reg
                                    + resize(to_signed(536870912, 32), 64);

                    -- FIX 4: gamma uses upd2_data + upd2_S_cur - newL
                    -- upd2_data = data_buf(t) - S_arr(si) from S_UPDATE
                    -- upd2_data + upd2_S_cur = data_buf(t)
                    -- data_buf(t) - newL = what gamma needs
                    dsp_gamma_x_dS <= gamma  * (upd2_data + upd2_S_cur
                                             - (dsp_alpha_x_data(61 downto 30)
                                             + dsp_1ma_x_LpT(61 downto 30)))
                                    + resize(to_signed(536870912, 32), 64);
                    dsp_1mg_x_Scur <= v_1mg  * upd2_S_cur
                                    + resize(to_signed(536870912, 32), 64);

                    -- Pass values through to S_UPDATE_WAIT
                    upd_fcast <= upd2_fcast;
                    upd_S_cur <= upd2_S_cur;
                    upd_si    <= upd2_si;
                    upd_t     <= upd2_t;
                    upd_data  <= upd2_data;

                    state <= S_UPDATE_WAIT;

                -- =============================================================
                -- S_UPDATE_WAIT: PIPELINE STAGE 4
                --
                -- Reads registered DSP results from S_UPDATE_MUL2.
                -- All values fully registered -- clean timing paths.
                --
                -- Path: read registered DSP + add + register = 5ns  OK
                -- =============================================================
                when S_UPDATE_WAIT =>

                    fitted_out   <= upd_fcast;
                    fitted_valid <= '1';

                    L_reg         <= upd_newL;
                    T_reg         <= dsp_beta_x_dL(61 downto 30)
                                   + dsp_1mb_x_T(61 downto 30);
                    S_arr(upd_si) <= dsp_gamma_x_dS(61 downto 30)
                                   + dsp_1mg_x_Scur(61 downto 30);

                    pipe_valid <= '0';

                    if upd_t = N_act - 1 then
                        -- si_counter was already incremented in S_UPDATE
                        -- for the last sample. It now points to the correct
                        -- starting seasonal index for forecasting.
                        -- Do NOT add 1 -- that causes fc_si_counter to start
                        -- at 1 instead of 0, shifting all seasonal factors
                        -- by one position and corrupting every forecast.
                        fc_si_counter <= si_counter;
                        k_cnt <= 0;
                        state <= S_FORECAST;
                    else
                        state <= S_UPDATE;
                    end if;

                -- =============================================================
                -- S_FORECAST: unchanged -- not on critical path
                -- =============================================================
                when S_FORECAST =>
                    k1_fp := to_signed((k_cnt + 1) * 65536, 32);

                    forecast_out   <= L_reg
                                    + fp_mul_data(k1_fp, T_reg)
                                    + S_arr(fc_si_counter);
                    forecast_valid <= '1';

                    if fc_si_counter = M_act - 1 then
                        fc_si_counter <= 0;
                    else
                        fc_si_counter <= fc_si_counter + 1;
                    end if;

                    if k_cnt = HORIZON_act - 1 then
                        last_forecast <= '1';
                        state         <= S_IDLE;
                    else
                        k_cnt <= k_cnt + 1;
                    end if;

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
