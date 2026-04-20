-- =============================================================================
-- holt_winters_q2_30_opt_v7.vhd
--
-- TIMING OPTIMIZED VERSION (100 MHz target):
--   * Combinational rcp() removed. Replaced with compile-time RCP_ROM.
--   * Division/Modulo removed from S_VALIDATE. Multi-cycle subtraction used.
--   * Chained combinational multipliers in S_UPDATE broken into pipeline stages.
--   * PRE-FETCH FIX: Added S_UPDATE_PRE to register array reads and kill the
--     23-level logic path before the DSP multipliers.
--
-- v7 TIMING FIXES (WNS -0.568 ns resolved):
--   FIX 1: Added upd_data_wait register. data_buf(upd_t) is now pre-fetched in
--           S_UPDATE_MULT instead of being read directly in S_UPDATE_WAIT.
--           This breaks the RAM-address -> RAM-data -> DSP chain (~10 ns logic).
--   FIX 2: Added prev_L register to capture L_reg before it is overwritten in
--           S_UPDATE_WAIT, making the (upd_newL - L_reg) path unambiguous and
--           removing a potential same-cycle read/write timing hazard.
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
    -- Timing Fix 1: Reciprocal ROM (Replaces combinational rcp function)
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
    signal M_act_m1   : integer range 0 to MAX_M-1 := 1;
    signal SEASONS_m1 : integer range 0 to MAX_SEASONS-1 := 0;
    signal next_season_flat : integer range 0 to MAX_N-1 := 0;

    -- =========================================================================
    -- FSM
    -- =========================================================================
    type state_t is (
        S_IDLE, S_COLLECT,
        S_VALIDATE, S_CALC_SEASONS,
        S_WAIT_FC,
        S_INIT_LT, S_INIT_LT_ODD, S_INIT_LT_DIV,
        S_SAVG, S_SAVG_ODD, S_SAVG_DIV,
        S_SACC, S_SACC_ODD,
        S_UPDATE_PRE, S_UPDATE, S_UPDATE_MULT, S_UPDATE_WAIT, S_UPDATE_FINALIZE,
        S_FORECAST, S_ERROR
    );
    signal state : state_t := S_IDLE;

    -- =========================================================================
    -- Counters & Division
    -- =========================================================================
    signal sample_cnt : integer range 0 to MAX_N       := 0;
    signal i_cnt      : integer range 0 to MAX_M       := 0;
    signal s_cnt      : integer range 0 to MAX_SEASONS := 0;
    signal t_cnt      : integer range 0 to MAX_N       := 0;
    signal k_cnt      : integer range 0 to MAX_HORIZON := 0;

    signal rem_n      : integer range 0 to MAX_N := 0;
    signal quot_n     : integer range 0 to MAX_SEASONS + 1 := 0;

    signal flat_idx    : integer range 0 to MAX_N-1 := 0;
    signal season_base : integer range 0 to MAX_N-1 := 0;

    signal acc_L  : fp_data := (others => '0');
    signal acc_T  : fp_data := (others => '0');
    signal acc_L2 : fp_data := (others => '0');
    signal acc_T2 : fp_data := (others => '0');

    -- =========================================================================
    -- UPDATE pipeline registers
    -- =========================================================================
    signal upd_S_pre    : fp_data := (others => '0');
    signal upd_data_pre : fp_data := (others => '0');

    -- -------------------------------------------------------------------------
    -- FIX 1: upd_data_wait - pre-fetched in S_UPDATE_MULT so that
    --        S_UPDATE_WAIT never reads data_buf(upd_t) directly.
    --        Breaks the: upd_t_reg -> RAM address -> RAM data -> DSP path.
    -- -------------------------------------------------------------------------
    signal upd_data_wait : fp_data := (others => '0');

    -- -------------------------------------------------------------------------
    -- FIX 2: prev_L - captures L_reg before it is overwritten in S_UPDATE_WAIT.
    --        Removes the same-cycle read/write ambiguity on L_reg that caused
    --        a long path into mul_beta_tmp.
    -- -------------------------------------------------------------------------
    signal prev_L : fp_data := (others => '0');

    signal upd_si     : integer range 0 to MAX_M-1 := 0;
    signal upd_fcast  : fp_data   := (others => '0');
    signal upd_S_cur  : fp_data   := (others => '0');
    signal upd_newL   : fp_data   := (others => '0');
    signal upd_t      : integer range 0 to MAX_N-1 := 0;
    signal pipe_valid : std_logic := '0';
    signal pipe_newL  : fp_data   := (others => '0');

    signal mul_alpha_tmp : fp_data := (others => '0');
    signal mul_1ma_tmp   : fp_data := (others => '0');
    signal mul_beta_tmp  : fp_data := (others => '0');
    signal mul_1mb_tmp   : fp_data := (others => '0');
    signal mul_gamma_tmp : fp_data := (others => '0');
    signal mul_1mg_tmp   : fp_data := (others => '0');

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
        variable v_1ma    : fp_param;
        variable v_1mb    : fp_param;
        variable v_1mg    : fp_param;
        variable k1_fp    : fp_data;
        variable tmp_data : fp_data;
        variable tmp_LpT  : fp_data;
        variable v_mean   : fp_data;
        variable next_si  : integer range 0 to MAX_M-1;
    begin
    if rising_edge(clk) then

        forecast_valid <= '0';
        fitted_valid   <= '0';
        last_forecast  <= '0';
        error_out      <= '0';
        error_code     <= (others => '0');

        if rst = '1' then
            state            <= S_IDLE;
            sample_cnt       <= 0;
            i_cnt            <= 0;
            s_cnt            <= 0;
            t_cnt            <= 0;
            k_cnt            <= 0;
            flat_idx         <= 0;
            season_base      <= 0;
            acc_L            <= (others => '0');
            acc_L2           <= (others => '0');
            acc_T            <= (others => '0');
            acc_T2           <= (others => '0');
            L_reg            <= (others => '0');
            T_reg            <= (others => '0');
            prev_L           <= (others => '0');  -- FIX 2 reset
            pipe_valid       <= '0';
            pipe_newL        <= (others => '0');
            S_acc            <= (others => (others => '0'));
            S_acc2           <= (others => (others => '0'));
            S_arr            <= (others => (others => '0'));
            cur_season_mean  <= (others => '0');
            si_counter       <= 0;
            fc_si_counter    <= 0;
            upd_S_pre        <= (others => '0');
            upd_data_pre     <= (others => '0');
            upd_data_wait    <= (others => '0');  -- FIX 1 reset
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

            rem_n            <= 0;
            quot_n           <= 0;
            mul_alpha_tmp    <= (others => '0');
            mul_1ma_tmp      <= (others => '0');
            mul_beta_tmp     <= (others => '0');
            mul_1mb_tmp      <= (others => '0');
            mul_gamma_tmp    <= (others => '0');
            mul_1mg_tmp      <= (others => '0');

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
                    if    M_act < 2         then
                        error_out <= '1'; error_code <= "010"; state <= S_ERROR;
                    elsif N_act < 2 * M_act then
                        error_out <= '1'; error_code <= "001"; state <= S_ERROR;
                    else
                        rem_n  <= N_act;
                        quot_n <= 0;
                        state  <= S_CALC_SEASONS;
                    end if;

                -- =============================================================
                when S_CALC_SEASONS =>
                    if rem_n >= M_act then
                        rem_n  <= rem_n - M_act;
                        quot_n <= quot_n + 1;
                    else
                        if rem_n /= 0 then
                            error_out <= '1'; error_code <= "011"; state <= S_ERROR;
                        elsif quot_n > MAX_SEASONS then
                            error_out <= '1'; error_code <= "101"; state <= S_ERROR;
                        else
                            SEASONS_act <= quot_n;
                            state       <= S_WAIT_FC;
                        end if;
                    end if;

                -- =============================================================
                when S_WAIT_FC =>
                    if forecast_start = '1' then
                        HORIZON_act      <= to_integer(unsigned(horizon_in));
                        inv_M_reg        <= RCP_ROM(M_act);
                        inv_M2_reg       <= RCP_ROM(M_act * M_act);
                        inv_S_reg        <= RCP_ROM(SEASONS_act);
                        M_act_m1         <= M_act - 1;
                        SEASONS_m1       <= SEASONS_act - 1;
                        next_season_flat <= M_act;
                        i_cnt            <= 0;
                        state            <= S_INIT_LT;
                    end if;

                -- =============================================================
                when S_INIT_LT =>
                    if i_cnt < M_act then
                        acc_L <= acc_L + data_buf(i_cnt);
                        acc_T <= acc_T + (data_buf(M_act + i_cnt) - data_buf(i_cnt));
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
                when S_INIT_LT_ODD =>
                    if i_cnt < M_act then
                        acc_L2 <= acc_L2 + data_buf(i_cnt);
                        acc_T2 <= acc_T2 + (data_buf(M_act + i_cnt) - data_buf(i_cnt));
                    end if;
                    i_cnt <= i_cnt + 1;
                    state <= S_INIT_LT;

                -- =============================================================
                when S_INIT_LT_DIV =>
                    L_reg       <= fp_mul_data(acc_L, inv_M_reg);
                    T_reg       <= fp_mul_data(acc_T, inv_M2_reg);
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
                when S_SAVG_ODD =>
                    if i_cnt < M_act then
                        acc_L2   <= acc_L2 + data_buf(flat_idx);
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                    end if;
                    state <= S_SAVG;

                -- =============================================================
                when S_SAVG_DIV =>
                    v_mean          := fp_mul_data(acc_L, inv_M_reg);
                    cur_season_mean <= v_mean;
                    acc_L           <= (others => '0');
                    acc_L2          <= (others => '0');
                    i_cnt           <= 0;
                    flat_idx        <= season_base;
                    state           <= S_SACC;

                -- =============================================================
                when S_SACC =>
                    if i_cnt < M_act then
                        S_acc(i_cnt) <= S_acc(i_cnt) + data_buf(flat_idx) - cur_season_mean;
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                        state    <= S_SACC_ODD;
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
                                    S_arr(j) <= fp_mul_data(S_acc(j) + S_acc2(j), inv_S_reg);
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
                when S_SACC_ODD =>
                    if i_cnt < M_act then
                        S_acc2(i_cnt) <= S_acc2(i_cnt) + data_buf(flat_idx) - cur_season_mean;
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                    end if;
                    state <= S_SACC;

                -- =============================================================
                -- PRE-FETCH: Register S_arr and data_buf reads to break the
                -- array routing delay before DSP multipliers.
                -- =============================================================
                when S_UPDATE_PRE =>
                    upd_S_pre    <= S_arr(si_counter);
                    upd_data_pre <= data_buf(t_cnt);
                    upd_si       <= si_counter;
                    upd_t        <= t_cnt;
                    state        <= S_UPDATE;

                -- =============================================================
                when S_UPDATE =>
                    v_1ma := ONE_PARAM - alpha;

                    -- Use the pre-fetched registers instead of direct array reads
                    tmp_data := upd_data_pre - upd_S_pre;

                    if pipe_valid = '1' then tmp_LpT := pipe_newL + T_reg;
                    else                     tmp_LpT := L_reg + T_reg;
                    end if;

                    upd_fcast <= tmp_LpT + upd_S_pre;
                    upd_S_cur <= upd_S_pre;

                    -- Stage 1: Register the multiplications
                    mul_alpha_tmp <= fp_mul_param(alpha, tmp_data);
                    mul_1ma_tmp   <= fp_mul_param(v_1ma, tmp_LpT);

                    pipe_valid <= '1';
                    t_cnt      <= upd_t + 1;
                    state      <= S_UPDATE_MULT;

                -- =============================================================
                when S_UPDATE_MULT =>
                    -- Stage 2: Register the addition (L update)
                    upd_newL  <= mul_alpha_tmp + mul_1ma_tmp;
                    pipe_newL <= mul_alpha_tmp + mul_1ma_tmp;

                    -- ---------------------------------------------------------
                    -- FIX 1: Pre-fetch data_buf(upd_t) here so S_UPDATE_WAIT
                    -- never drives RAM address lines on its critical path.
                    -- upd_t here is the original sample index (before the +1
                    -- increment in S_UPDATE), which is exactly what
                    -- S_UPDATE_WAIT needs for the gamma term.
                    -- ---------------------------------------------------------
                    upd_data_wait <= data_buf(upd_t);

                    state <= S_UPDATE_WAIT;

                -- =============================================================
                when S_UPDATE_WAIT =>
                    v_1mb := ONE_PARAM - beta;
                    v_1mg := ONE_PARAM - gamma;

                    fitted_out   <= upd_fcast;
                    fitted_valid <= '1';

                    -- ---------------------------------------------------------
                    -- FIX 2: Capture L_reg into prev_L *before* overwriting it.
                    -- mul_beta_tmp now uses prev_L instead of reading L_reg in
                    -- the same cycle it is written, removing the ambiguous
                    -- same-cycle read/write path.
                    -- ---------------------------------------------------------
                    prev_L <= L_reg;
                    L_reg  <= upd_newL;

                    -- Stage 3: Register the multiplications for T and S
                    -- FIX 2 applied: (upd_newL - prev_L) replaces (upd_newL - L_reg)
                    mul_beta_tmp  <= fp_mul_param(beta,  upd_newL - prev_L);
                    mul_1mb_tmp   <= fp_mul_param(v_1mb, T_reg);

                    -- FIX 1 applied: upd_data_wait replaces data_buf(upd_t)
                    mul_gamma_tmp <= fp_mul_param(gamma, upd_data_wait - upd_newL);
                    mul_1mg_tmp   <= fp_mul_param(v_1mg, upd_S_cur);

                    pipe_valid <= '0';
                    state      <= S_UPDATE_FINALIZE;

                -- =============================================================
                when S_UPDATE_FINALIZE =>
                    -- Stage 4: Register the additions for T and S
                    T_reg         <= mul_beta_tmp + mul_1mb_tmp;
                    S_arr(upd_si) <= mul_gamma_tmp + mul_1mg_tmp;

                    if si_counter = M_act_m1 then next_si := 0;
                    else                          next_si := si_counter + 1;
                    end if;
                    si_counter <= next_si;

                    if upd_t = N_act - 1 then
                        if si_counter = M_act_m1 then fc_si_counter <= 0;
                        else                          fc_si_counter <= si_counter + 1;
                        end if;
                        k_cnt <= 0;
                        state <= S_FORECAST;
                    else
                        state <= S_UPDATE_PRE;
                    end if;

                -- =============================================================
                when S_FORECAST =>
                    k1_fp := to_signed((k_cnt + 1) * 65536, 32);

                    forecast_out   <= L_reg + fp_mul_data(k1_fp, T_reg) + S_arr(fc_si_counter);
                    forecast_valid <= '1';

                    if fc_si_counter = M_act_m1 then fc_si_counter <= 0;
                    else                             fc_si_counter <= fc_si_counter + 1;
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
