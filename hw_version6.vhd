-- =============================================================================
-- holt_winters_q2_30_opt_v4.vhd
--
-- ALL CHANGES FROM VERSION 3 RETAINED PLUS:
--
-- FIX: S_SAVG / S_SACC single address expression via flat_idx
--
--   PROBLEM (Version 3):
--     S_SAVG and S_SACC used two address expressions:
--       data_buf(idx)     where idx = s_cnt * M_act + i_cnt  -- expression 1
--       data_buf(idx + 1)                                     -- expression 2
--     The +1 offset (p_4_in nets) forced an extra RAM replica.
--
--   SOLUTION:
--     Introduce flat_idx signal - a flat running pointer into data_buf.
--     Introduce season_base signal - tracks start of current season.
--     Split S_SAVG into S_SAVG / S_SAVG_ODD alternating states.
--     Split S_SACC into S_SACC / S_SACC_ODD alternating states.
--     Both even and odd states read data_buf(flat_idx) - ONE expression.
--     flat_idx incremented each cycle so no +1 offset ever appears.
--     season_base used for clean rewind in S_SAVG_DIV instead of
--     flat_idx - M_act subtraction which caused range errors.
--
--   BUG FIXED vs first V4 attempt:
--     flat_idx range 0 to MAX_N-1 caused range violation on rewind.
--     flat_idx - M_act rewind was unreliable for multi-season case.
--     Fix: season_base tracks current season start explicitly.
--          S_SAVG_DIV rewinds flat_idx to season_base (always valid).
--          S_SACC advances season_base by M_act on season transition.
--
--   ADDRESS EXPRESSIONS ACROSS ALL STATES:
--     S_INIT_LT / S_INIT_LT_ODD : data_buf(i_cnt), data_buf(M_act+i_cnt)
--     S_SAVG / S_SAVG_ODD        : data_buf(flat_idx)
--     S_SACC / S_SACC_ODD        : data_buf(flat_idx)
--     S_UPDATE                   : data_buf(t_cnt)
--     S_UPDATE_WAIT              : data_buf(upd_t)
--
--   UNIQUE ADDRESS EXPRESSIONS: 4
--     1. i_cnt
--     2. M_act + i_cnt
--     3. flat_idx  (shared across SAVG, SAVG_ODD, SACC, SACC_ODD)
--     4. t_cnt
--     5. upd_t
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

    signal data_buf : data_array := (others => (others => '0'));
    attribute ram_style : string;
    attribute ram_style of data_buf : signal is "distributed";

    signal S_arr           : season_array := (others => (others => '0'));
    signal S_acc           : season_array := (others => (others => '0'));
    signal S_acc2          : season_array := (others => (others => '0'));
    signal cur_season_mean : fp_data      := (others => '0');

    signal L_reg : fp_data := (others => '0');
    signal T_reg : fp_data := (others => '0');

    -- =========================================================================
    -- FSM states
    -- =========================================================================
    type state_t is (
        S_IDLE, S_COLLECT, S_VALIDATE, S_WAIT_FC,
        S_INIT_LT,
        S_INIT_LT_ODD,
        S_INIT_LT_DIV,
        S_SAVG,
        S_SAVG_ODD,
        S_SAVG_DIV,
        S_SACC,
        S_SACC_ODD,
        S_UPDATE, S_UPDATE_WAIT,
        S_FORECAST, S_ERROR
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
    -- flat_idx: flat running pointer into data_buf for SAVG and SACC phases
    --
    -- Replaces idx = s_cnt * M_act + i_cnt and idx+1 expressions.
    -- data_buf(flat_idx) is the single address expression used across
    -- S_SAVG, S_SAVG_ODD, S_SACC, and S_SACC_ODD.
    --
    -- season_base: tracks the start index of the current season.
    -- Used by S_SAVG_DIV to rewind flat_idx back to season start
    -- so S_SACC re-reads the same data that S_SAVG just averaged.
    -- Avoids flat_idx - M_act subtraction which caused range violations.
    -- =========================================================================
    signal flat_idx    : integer range 0 to MAX_N-1 := 0;
    signal season_base : integer range 0 to MAX_N-1 := 0;

    -- =========================================================================
    -- Primary accumulators
    -- =========================================================================
    signal acc_L : fp_data := (others => '0');
    signal acc_T : fp_data := (others => '0');

    -- =========================================================================
    -- Secondary accumulators
    -- =========================================================================
    signal acc_L2 : fp_data := (others => '0');
    signal acc_T2 : fp_data := (others => '0');

    -- =========================================================================
    -- UPDATE pipeline registers
    -- =========================================================================
    signal upd_si     : integer range 0 to MAX_M-1 := 0;
    signal upd_fcast  : fp_data   := (others => '0');
    signal upd_S_cur  : fp_data   := (others => '0');
    signal upd_newL   : fp_data   := (others => '0');
    signal upd_t      : integer range 0 to MAX_N-1 := 0;
    signal pipe_valid : std_logic := '0';
    signal pipe_newL  : fp_data   := (others => '0');

    -- =========================================================================
    -- Season index counters
    -- =========================================================================
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

    -- =========================================================================
    -- rcp: 1/x in Q16.16
    -- =========================================================================
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
        variable k1_fp       : fp_data;
        variable tmp_data    : fp_data;
        variable tmp_LpT     : fp_data;
        variable last_sacc   : fp_data;
        variable last_sacc2  : fp_data;
        variable v_mean      : fp_data;
        variable next_si     : integer range 0 to MAX_M-1;
    begin
    if rising_edge(clk) then

        forecast_valid <= '0';
        fitted_valid   <= '0';
        last_forecast  <= '0';
        error_out      <= '0';
        error_code     <= (others => '0');

        if rst = '1' then
            state           <= S_IDLE;
            sample_cnt      <= 0;
            i_cnt           <= 0;
            s_cnt           <= 0;
            t_cnt           <= 0;
            k_cnt           <= 0;
            flat_idx        <= 0;
            season_base     <= 0;
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
            upd_si          <= 0;
            upd_fcast       <= (others => '0');
            upd_S_cur       <= (others => '0');
            upd_newL        <= (others => '0');
            upd_t           <= 0;

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
                            N_act      <= MAX_N;
                            sample_cnt <= 0;
                            acc_L      <= (others => '0');
                            acc_L2     <= (others => '0');
                            acc_T      <= (others => '0');
                            acc_T2     <= (others => '0');
                            i_cnt      <= 0;
                            state      <= S_VALIDATE;
                        else
                            sample_cnt <= sample_cnt + 1;
                        end if;
                    elsif sample_cnt > 0 then
                        N_act      <= sample_cnt;
                        sample_cnt <= 0;
                        acc_L      <= (others => '0');
                        acc_L2     <= (others => '0');
                        acc_T      <= (others => '0');
                        acc_T2     <= (others => '0');
                        i_cnt      <= 0;
                        state      <= S_VALIDATE;
                    end if;

                -- =============================================================
                when S_VALIDATE =>
                    if    M_act < 2                     then
                        error_out <= '1'; error_code <= "010";
                        state     <= S_ERROR;
                    elsif N_act < 2 * M_act             then
                        error_out <= '1'; error_code <= "001";
                        state     <= S_ERROR;
                    elsif not divides(M_act, N_act)     then
                        error_out <= '1'; error_code <= "011";
                        state     <= S_ERROR;
                    elsif (N_act / M_act) > MAX_SEASONS then
                        error_out <= '1'; error_code <= "101";
                        state     <= S_ERROR;
                    else
                        SEASONS_act <= N_act / M_act;
                        state       <= S_WAIT_FC;
                    end if;

                -- =============================================================
                when S_WAIT_FC =>
                    if forecast_start = '1' then
                        HORIZON_act <= to_integer(unsigned(horizon_in));
                        i_cnt       <= 0;
                        state       <= S_INIT_LT;
                    end if;

                -- =============================================================
                -- S_INIT_LT: reads data_buf(i_cnt) and data_buf(M_act+i_cnt)
                --
                -- i_cnt incremented here before S_INIT_LT_ODD so that
                -- S_INIT_LT_ODD reuses the same two address expressions.
                -- =============================================================
                when S_INIT_LT =>
                    if i_cnt < M_act then
                        acc_L <= acc_L + data_buf(i_cnt);
                        acc_T <= acc_T + (data_buf(M_act + i_cnt)
                                        - data_buf(i_cnt));
                        i_cnt <= i_cnt + 1;
                        state <= S_INIT_LT_ODD;
                    else
                        acc_L  <= acc_L + acc_L2;
                        acc_T  <= acc_T + acc_T2;
                        acc_L2 <= (others => '0');
                        acc_T2 <= (others => '0');
                        i_cnt  <= 0;
                        state  <= S_INIT_LT_DIV;
                    end if;

                -- =============================================================
                -- S_INIT_LT_ODD: same expressions as S_INIT_LT
                -- i_cnt already incremented so points to odd element
                -- =============================================================
                when S_INIT_LT_ODD =>
                    if i_cnt < M_act then
                        acc_L2 <= acc_L2 + data_buf(i_cnt);
                        acc_T2 <= acc_T2 + (data_buf(M_act + i_cnt)
                                          - data_buf(i_cnt));
                    end if;
                    i_cnt <= i_cnt + 1;
                    state <= S_INIT_LT;

                -- =============================================================
                -- S_INIT_LT_DIV: compute L0 and T0, no data_buf reads
                --
                -- flat_idx and season_base both set to 0 here.
                -- Season 0 starts at data_buf(0).
                -- =============================================================
                when S_INIT_LT_DIV =>
                    inv_M       := rcp(M_act);
                    inv_M2      := rcp(M_act * M_act);
                    L_reg       <= fp_mul_data(acc_L, inv_M);
                    T_reg       <= fp_mul_data(acc_T, inv_M2);
                    S_acc       <= (others => (others => '0'));
                    S_acc2      <= (others => (others => '0'));
                    acc_L       <= (others => '0');
                    acc_L2      <= (others => '0');
                    acc_T       <= (others => '0');
                    acc_T2      <= (others => '0');
                    i_cnt       <= 0;
                    s_cnt       <= 0;
                    flat_idx    <= 0;
                    season_base <= 0;
                    state       <= S_SAVG;

                -- =============================================================
                -- S_SAVG: season mean accumulation - even element
                --
                -- Reads data_buf(flat_idx) - single address expression.
                -- Increments flat_idx and i_cnt then goes to S_SAVG_ODD.
                -- When i_cnt reaches M_act all elements read - go to DIV.
                -- =============================================================
                when S_SAVG =>
                    if i_cnt < M_act then
                        acc_L    <= acc_L + data_buf(flat_idx);
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                        state    <= S_SAVG_ODD;
                    else
                        acc_L  <= acc_L + acc_L2;
                        acc_L2 <= (others => '0');
                        i_cnt  <= 0;
                        state  <= S_SAVG_DIV;
                    end if;

                -- =============================================================
                -- S_SAVG_ODD: season mean accumulation - odd element
                --
                -- Reads data_buf(flat_idx) - SAME expression as S_SAVG.
                -- flat_idx was already incremented in S_SAVG so it now
                -- points to the odd element. No +1 offset needed.
                -- Skips read if i_cnt >= M_act (odd M_act case).
                -- =============================================================
                when S_SAVG_ODD =>
                    if i_cnt < M_act then
                        acc_L2   <= acc_L2 + data_buf(flat_idx);
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                    end if;
                    state <= S_SAVG;

                -- =============================================================
                -- S_SAVG_DIV: compute season mean, no data_buf reads
                --
                -- flat_idx rewound to season_base so S_SACC re-reads
                -- the same season elements that S_SAVG just averaged.
                -- Using season_base avoids range violation from subtraction.
                -- =============================================================
                when S_SAVG_DIV =>
                    inv_M           := rcp(M_act);
                    v_mean          := fp_mul_data(acc_L, inv_M);
                    cur_season_mean <= v_mean;
                    acc_L           <= (others => '0');
                    acc_L2          <= (others => '0');
                    i_cnt           <= 0;
                    flat_idx        <= season_base;
                    state           <= S_SACC;

                -- =============================================================
                -- S_SACC: seasonal deviation accumulation - even element
                --
                -- Reads data_buf(flat_idx) - SAME expression as S_SAVG.
                -- Increments flat_idx and i_cnt then goes to S_SACC_ODD.
                --
                -- When i_cnt reaches M_act season is complete:
                --   More seasons: advance season_base by M_act,
                --                 set flat_idx to new season_base,
                --                 reset i_cnt, go to S_SAVG
                --   Last season:  write S_arr and go to S_UPDATE
                -- =============================================================
                when S_SACC =>
                    if i_cnt < M_act then
                        S_acc(i_cnt) <= S_acc(i_cnt)
                                      + data_buf(flat_idx)
                                      - cur_season_mean;
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                        state    <= S_SACC_ODD;
                    else
                        if s_cnt < SEASONS_act - 1 then
                            -- Advance to next season
                            -- season_base advances by M_act
                            -- flat_idx set to new season_base for S_SAVG
                            s_cnt       <= s_cnt + 1;
                            i_cnt       <= 0;
                            season_base <= season_base + M_act;
                            flat_idx    <= season_base + M_act;
                            state       <= S_SAVG;
                        else
                            -- All seasons done
                            -- Combine S_acc and S_acc2 and write S_arr
                            inv_SEASONS := rcp(SEASONS_act);
                            for j in 0 to MAX_M - 1 loop
                                if j < M_act then
                                    S_arr(j) <= fp_mul_data(
                                        S_acc(j) + S_acc2(j),
                                        inv_SEASONS);
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
                -- S_SACC_ODD: seasonal deviation accumulation - odd element
                --
                -- Reads data_buf(flat_idx) - SAME expression as S_SACC.
                -- flat_idx was already incremented in S_SACC.
                -- Skips read if i_cnt >= M_act (odd M_act case).
                -- =============================================================
                when S_SACC_ODD =>
                    if i_cnt < M_act then
                        S_acc2(i_cnt) <= S_acc2(i_cnt)
                                       + data_buf(flat_idx)
                                       - cur_season_mean;
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                    end if;
                    state <= S_SACC;

                -- =============================================================
                -- S_UPDATE: pipeline stage 1
                -- Single read: data_buf(t_cnt)
                -- =============================================================
                when S_UPDATE =>
                    v_si  := si_counter;
                    v_1ma := ONE_PARAM - alpha;

                    tmp_data := data_buf(t_cnt) - S_arr(v_si);

                    if pipe_valid = '1' then
                        tmp_LpT := pipe_newL + T_reg;
                    else
                        tmp_LpT := L_reg + T_reg;
                    end if;

                    upd_fcast <= tmp_LpT + S_arr(v_si);
                    upd_newL  <= fp_mul_param(alpha, tmp_data)
                               + fp_mul_param(v_1ma, tmp_LpT);
                    upd_S_cur <= S_arr(v_si);
                    upd_si    <= v_si;
                    upd_t     <= t_cnt;

                    if pipe_valid = '1' then
                        pipe_newL <= fp_mul_param(alpha, tmp_data)
                                   + fp_mul_param(v_1ma, pipe_newL + T_reg);
                    else
                        pipe_newL <= fp_mul_param(alpha, tmp_data)
                                   + fp_mul_param(v_1ma, L_reg + T_reg);
                    end if;

                    pipe_valid <= '1';
                    t_cnt      <= t_cnt + 1;
                    state      <= S_UPDATE_WAIT;

                -- =============================================================
                -- S_UPDATE_WAIT: pipeline stage 2
                -- =============================================================
                when S_UPDATE_WAIT =>
                    v_1mb := ONE_PARAM - beta;
                    v_1mg := ONE_PARAM - gamma;

                    fitted_out   <= upd_fcast;
                    fitted_valid <= '1';

                    L_reg         <= upd_newL;
                    T_reg         <= fp_mul_param(beta,  upd_newL - L_reg)
                                   + fp_mul_param(v_1mb, T_reg);
                    S_arr(upd_si) <= fp_mul_param(gamma,
                                        data_buf(upd_t) - upd_newL)
                                   + fp_mul_param(v_1mg, upd_S_cur);

                    pipe_valid <= '0';

                    if si_counter = M_act - 1 then
                        next_si := 0;
                    else
                        next_si := si_counter + 1;
                    end if;
                    si_counter <= next_si;

                    if upd_t = N_act - 1 then
                        if si_counter = M_act - 1 then
                            fc_si_counter <= 0;
                        else
                            fc_si_counter <= si_counter + 1;
                        end if;
                        k_cnt <= 0;
                        state <= S_FORECAST;
                    else
                        state <= S_UPDATE;
                    end if;

                -- =============================================================
                -- S_FORECAST: no data_buf reads
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
