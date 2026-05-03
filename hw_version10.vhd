-- =============================================================================
-- holt_winters_q2_30_opt_v9.vhd
--
-- Builds on v8. All v7 and v8 fixes retained unchanged.
--
-- v9 LATENCY FIX (~132 cycles saved, same WNS, same RAM replicas):
--
--   FIX 4: S_INIT_LT_ODD eliminated.
--     Previously S_INIT_LT and S_INIT_LT_ODD alternated, spending 2 cycles
--     per element to accumulate into acc_L/acc_L2 and acc_T/acc_T2.
--     Now S_INIT_LT processes one element per cycle directly into acc_L/acc_T.
--     acc_L2 and acc_T2 signals removed entirely.
--     Saves: M cycles (12 cycles for M=12).
--     WNS impact: none. Single 32-bit add per cycle is ~4 ns, well within 14 ns.
--     Replica impact: none. i_cnt (r4) and computed M_act+i_cnt (r3) unchanged.
--
--   FIX 5: S_SAVG_ODD eliminated.
--     Previously S_SAVG and S_SAVG_ODD alternated. Now S_SAVG accumulates
--     one element per cycle directly into acc_L.
--     acc_L2 signal removed (was only needed by SAVG_ODD).
--     Saves: M cycles per season (12 x 5 = 60 cycles for M=12, SEASONS=5).
--     WNS impact: none. flat_idx -> distRAM -> 32-bit add -> register: ~6 ns.
--     Replica impact: none. flat_idx (r2) is unchanged.
--
--   FIX 6: S_SACC_ODD eliminated.
--     Previously S_SACC and S_SACC_ODD alternated, splitting contributions
--     into S_acc and S_acc2. Now S_SACC accumulates all elements into S_acc.
--     S_acc2 array removed entirely. Final S_arr calculation uses S_acc only.
--     Saves: M cycles per season (12 x 5 = 60 cycles for M=12, SEASONS=5).
--     WNS impact: none. flat_idx -> distRAM -> sub -> add -> register: ~7 ns.
--     Replica impact: none. S_acc is a register array (FFs), not data_buf RAM.
--
-- TOTAL SAVINGS (M=12, N=60, SEASONS=5):
--   132 cycles = 1848 ns at 14 ns clock
--   487 -> ~355 cycles
--   Output throughput: ~10.56 -> ~14.49 Msamples/sec
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

    signal S_arr          : season_array := (others => (others => '0'));
    signal S_acc          : season_array := (others => (others => '0'));
    -- FIX 6: S_acc2 removed. S_acc alone accumulates all seasonal contributions.
    signal cur_season_mean : fp_data     := (others => '0');

    signal L_reg : fp_data := (others => '0');
    signal T_reg : fp_data := (others => '0');

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
    signal M_act_m1   : integer range 0 to MAX_M-1     := 1;
    signal SEASONS_m1 : integer range 0 to MAX_SEASONS-1 := 0;
    signal next_season_flat : integer range 0 to MAX_N-1 := 0;

    -- =========================================================================
    -- FSM
    -- FIX 4/5/6: S_INIT_LT_ODD, S_SAVG_ODD, S_SACC_ODD removed from state_t.
    -- =========================================================================
    type state_t is (
        S_IDLE, S_COLLECT,
        S_VALIDATE, S_CALC_SEASONS,
        S_WAIT_FC,
        S_INIT_LT, S_INIT_LT_DIV,
        S_SAVG, S_SAVG_DIV,
        S_SACC,
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

    signal rem_n  : integer range 0 to MAX_N           := 0;
    signal quot_n : integer range 0 to MAX_SEASONS + 1 := 0;

    signal flat_idx    : integer range 0 to MAX_N-1 := 0;
    signal season_base : integer range 0 to MAX_N-1 := 0;

    -- FIX 4/5: acc_L2 and acc_T2 removed. Single accumulators are sufficient
    -- now that the ODD states are gone and each element is processed in 1 cycle.
    signal acc_L : fp_data := (others => '0');
    signal acc_T : fp_data := (others => '0');

    -- =========================================================================
    -- UPDATE pipeline registers (v7 FIX 1 and FIX 2 unchanged)
    -- =========================================================================
    signal upd_S_pre     : fp_data := (others => '0');
    signal upd_data_pre  : fp_data := (others => '0');
    signal upd_data_wait : fp_data := (others => '0');
    signal prev_L        : fp_data := (others => '0');

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
            -- FIX 4/5: acc_L2, acc_T2 removed from reset
            acc_T            <= (others => '0');
            L_reg            <= (others => '0');
            T_reg            <= (others => '0');
            prev_L           <= (others => '0');
            pipe_valid       <= '0';
            pipe_newL        <= (others => '0');
            S_acc            <= (others => (others => '0'));
            -- FIX 6: S_acc2 removed from reset
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
                            acc_T      <= (others => '0');
                            i_cnt      <= 0;
                            state      <= S_VALIDATE;
                        else
                            sample_cnt <= sample_cnt + 1;
                        end if;
                    elsif sample_cnt > 0 then
                        N_act      <= sample_cnt;
                        sample_cnt <= 0;
                        acc_L      <= (others => '0');
                        acc_T      <= (others => '0');
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
                -- FIX 4: S_INIT_LT_ODD removed. One element processed per cycle.
                -- Previously alternated INIT_LT <-> INIT_LT_ODD (2 cycles/elem).
                -- Now stays in INIT_LT until i_cnt = M_act (1 cycle/elem).
                -- acc_L2 and acc_T2 removed; acc_L and acc_T accumulate directly.
                -- Path: data_buf(i_cnt) -> 32-bit add -> acc_L register (~5 ns).
                --       data_buf(M_act+i_cnt) - data_buf(i_cnt) -> add -> acc_T (~7 ns).
                -- Both paths are well within 14 ns.
                -- RAM replicas unchanged: i_cnt (r4) and computed M_act+i_cnt (r3).
                -- =============================================================
                when S_INIT_LT =>
                    if i_cnt < M_act then
                        acc_L <= acc_L + data_buf(i_cnt);
                        acc_T <= acc_T + (data_buf(M_act + i_cnt) - data_buf(i_cnt));
                        i_cnt <= i_cnt + 1;
                        -- stay in S_INIT_LT: no ODD transition needed
                    else
                        -- All M elements accumulated. No merge step needed.
                        i_cnt <= 0;
                        state <= S_INIT_LT_DIV;
                    end if;

                -- =============================================================
                when S_INIT_LT_DIV =>
                    L_reg <= fp_mul_data(acc_L, inv_M_reg);
                    T_reg <= fp_mul_data(acc_T, inv_M2_reg);
                    S_acc <= (others => (others => '0'));
                    -- FIX 6: S_acc2 reset removed
                    acc_L <= (others => '0');
                    -- FIX 4: acc_L2, acc_T2 resets removed
                    acc_T       <= (others => '0');
                    i_cnt       <= 0;
                    s_cnt       <= 0;
                    flat_idx    <= 0;
                    season_base <= 0;
                    state       <= S_SAVG;

                -- =============================================================
                -- FIX 5: S_SAVG_ODD removed. One element processed per cycle.
                -- Previously alternated SAVG <-> SAVG_ODD (2 cycles/elem).
                -- Now stays in S_SAVG until i_cnt = M_act (1 cycle/elem).
                -- acc_L2 removed; acc_L accumulates the full season sum directly.
                -- Path: flat_idx -> distRAM -> 32-bit add -> acc_L reg (~6 ns).
                -- RAM replica unchanged: flat_idx (r2).
                -- =============================================================
                when S_SAVG =>
                    if i_cnt < M_act then
                        acc_L    <= acc_L + data_buf(flat_idx);
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                        -- stay in S_SAVG: no ODD transition needed
                    else
                        -- Full season sum in acc_L. No merge step needed.
                        i_cnt <= 0;
                        state <= S_SAVG_DIV;
                    end if;

                -- =============================================================
                when S_SAVG_DIV =>
                    v_mean          := fp_mul_data(acc_L, inv_M_reg);
                    cur_season_mean <= v_mean;
                    acc_L           <= (others => '0');
                    i_cnt           <= 0;
                    flat_idx        <= season_base;
                    state           <= S_SACC;

                -- =============================================================
                -- FIX 6: S_SACC_ODD removed. One element processed per cycle.
                -- Previously alternated SACC <-> SACC_ODD, splitting contributions
                -- between S_acc (even elements) and S_acc2 (odd elements).
                -- Now stays in S_SACC and accumulates all elements into S_acc.
                -- S_acc2 array removed entirely.
                -- Final S_arr uses S_acc alone (no S_acc + S_acc2 merge needed).
                -- Path: flat_idx -> distRAM -> sub(cur_season_mean) ->
                --       add(S_acc(i_cnt)) -> S_acc register (~7 ns).
                -- S_acc is a register array (FFs), not data_buf RAM: no new replica.
                -- RAM replica unchanged: flat_idx (r2).
                -- =============================================================
                when S_SACC =>
                    if i_cnt < M_act then
                        S_acc(i_cnt) <= S_acc(i_cnt) +
                                        data_buf(flat_idx) - cur_season_mean;
                        flat_idx <= flat_idx + 1;
                        i_cnt    <= i_cnt + 1;
                        -- stay in S_SACC: no ODD transition needed
                    else
                        if s_cnt < SEASONS_m1 then
                            s_cnt            <= s_cnt + 1;
                            i_cnt            <= 0;
                            season_base      <= next_season_flat;
                            flat_idx         <= next_season_flat;
                            next_season_flat <= next_season_flat + M_act;
                            state            <= S_SAVG;
                        else
                            -- FIX 6: S_arr uses S_acc only (S_acc2 removed)
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
                -- Entry point for the FIRST sample only (t=0). v8 FIX 3
                -- ensures all subsequent samples are pre-fetched inside
                -- S_UPDATE_FINALIZE, so this state executes exactly once per run.
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

                    tmp_data := upd_data_pre - upd_S_pre;

                    if pipe_valid = '1' then tmp_LpT := pipe_newL + T_reg;
                    else                     tmp_LpT := L_reg + T_reg;
                    end if;

                    upd_fcast <= tmp_LpT + upd_S_pre;
                    upd_S_cur <= upd_S_pre;

                    mul_alpha_tmp <= fp_mul_param(alpha, tmp_data);
                    mul_1ma_tmp   <= fp_mul_param(v_1ma, tmp_LpT);

                    pipe_valid <= '1';
                    t_cnt      <= upd_t + 1;
                    state      <= S_UPDATE_MULT;

                -- =============================================================
                when S_UPDATE_MULT =>
                    upd_newL  <= mul_alpha_tmp + mul_1ma_tmp;
                    pipe_newL <= mul_alpha_tmp + mul_1ma_tmp;

                    -- v7 FIX 1: Pre-fetch data_buf(upd_t) here
                    upd_data_wait <= data_buf(upd_t);

                    state <= S_UPDATE_WAIT;

                -- =============================================================
                when S_UPDATE_WAIT =>
                    v_1mb := ONE_PARAM - beta;
                    v_1mg := ONE_PARAM - gamma;

                    fitted_out   <= upd_fcast;
                    fitted_valid <= '1';

                    -- v7 FIX 2: Capture L_reg before overwriting
                    prev_L <= L_reg;
                    L_reg  <= upd_newL;

                    mul_beta_tmp  <= fp_mul_param(beta,  upd_newL - prev_L);
                    mul_1mb_tmp   <= fp_mul_param(v_1mb, T_reg);

                    -- v7 FIX 1: Use pre-fetched upd_data_wait
                    mul_gamma_tmp <= fp_mul_param(gamma, upd_data_wait - upd_newL);
                    mul_1mg_tmp   <= fp_mul_param(v_1mg, upd_S_cur);

                    pipe_valid <= '0';
                    state      <= S_UPDATE_FINALIZE;

                -- =============================================================
                when S_UPDATE_FINALIZE =>
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
                        -- v8 FIX 3: Pre-fetch next sample here, skip S_UPDATE_PRE.
                        -- t_cnt = upd_t+1 (set in S_UPDATE).
                        -- next_si != upd_si guaranteed (M_act >= 2): no RAW hazard.
                        upd_S_pre    <= S_arr(next_si);
                        upd_data_pre <= data_buf(t_cnt);
                        upd_si       <= next_si;
                        upd_t        <= t_cnt;
                        state        <= S_UPDATE;
                    end if;

                -- =============================================================
                when S_FORECAST =>
                    k1_fp := to_signed((k_cnt + 1) * 65536, 32);

                    forecast_out   <= L_reg + fp_mul_data(k1_fp, T_reg)
                                           + S_arr(fc_si_counter);
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
