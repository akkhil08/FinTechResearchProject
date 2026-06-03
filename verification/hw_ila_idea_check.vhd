-- =============================================================================
-- holt_winters.vhd  -  Holt-Winters Additive Exponential Smoothing
-- Stream-based, fully synthesisable for Xilinx ZedBoard (XC7Z020)
--
-- ALGORITHM (Additive Holt-Winters):
--   Level:   L(t) = alpha*(y(t)-S(t-m))   + (1-alpha)*(L(t-1)+T(t-1))
--   Trend:   T(t) = beta *(L(t)-L(t-1))   + (1-beta) *T(t-1)
--   Season:  S(t) = gamma*(y(t)-L(t))     + (1-gamma)*S(t-m)
--   Forecast:F(t+k)= L(t) + k*T(t) + S(t+k-m)
--
-- FIXED-POINT FORMAT: Q16.16 (signed 32-bit)
--   Real value = integer / 65536
--   ONE = 65536  (represents 1.0)
--
-- FSM STATES:
--   S_IDLE        -> wait for start pulse, latch M
--   S_COLLECT     -> stream data_in samples (valid_in counts N)
--   S_VALIDATE    -> check N/M constraints
--   S_WAIT_FC     -> wait for forecast_start to latch horizon
--   S_INIT_L      -> compute initial level L0 = mean(y[0..M-1])
--   S_INIT_T      -> compute initial trend T0
--   S_INIT_SAVG   -> compute seasonal averages per season
--   S_INIT_S      -> compute initial seasonal indices S[0..M-1]
--   S_UPDATE      -> recursive Holt-Winters update (stage 1)
--   S_UPDATE_WAIT -> recursive Holt-Winters update (stage 2, commit)
--   S_FORECAST    -> output HORIZON forecasts
--   S_ERROR       -> latch error, return to idle
--
-- ERROR CODES:
--   001 = N < 2*M         (need at least 2 complete seasons)
--   010 = M < 2           (season length must be >= 2)
--   011 = M does not divide N evenly
--   101 = N/M > MAX_SEASONS
-- =============================================================================

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity holt_winters is
    generic (
        MAX_N       : integer := 72;    -- maximum number of input samples
        MAX_M       : integer := 24;    -- maximum season length
        MAX_HORIZON : integer := 24;    -- maximum forecast horizon
        MAX_SEASONS : integer := 8      -- maximum number of complete seasons
    );
    port (
        clk             : in  std_logic;
        rst             : in  std_logic;

        -- Season length (latched on start pulse)
        m_in            : in  integer range 2 to MAX_M;

        -- Forecast horizon (latched on forecast_start pulse)
        forecast_start  : in  std_logic;
        horizon_in      : in  integer range 1 to MAX_HORIZON;

        -- Control
        start           : in  std_logic;
        alpha           : in  signed(31 downto 0);   -- Q16.16  e.g. 62259 = 0.95
        beta            : in  signed(31 downto 0);   -- Q16.16
        gamma           : in  signed(31 downto 0);   -- Q16.16

        -- Input data stream (Q16.16)
        data_in         : in  signed(31 downto 0);
        valid_in        : in  std_logic;

        -- Fitted output stream
        fitted_out      : out signed(31 downto 0);
        fitted_valid    : out std_logic;

        -- Forecast output stream
        forecast_out    : out signed(31 downto 0);
        forecast_valid  : out std_logic;
        last_forecast   : out std_logic;   -- pulses with final forecast_valid

        -- Error reporting
        error_out       : out std_logic;
        error_code      : out std_logic_vector(2 downto 0)
    );
end entity holt_winters;

architecture rtl of holt_winters is

    -- -------------------------------------------------------------------------
    -- Fixed-point type and constants
    -- -------------------------------------------------------------------------
    subtype fp is signed(31 downto 0);
    constant ONE : fp := to_signed(65536, 32);   -- 1.0 in Q16.16

    -- -------------------------------------------------------------------------
    -- Runtime parameters
    -- -------------------------------------------------------------------------
    signal N_act       : integer range 1 to MAX_N       := 1;
    signal M_act       : integer range 2 to MAX_M       := 2;
    signal HORIZON_act : integer range 1 to MAX_HORIZON := 1;
    signal SEASONS_act : integer range 1 to MAX_SEASONS := 1;

    -- -------------------------------------------------------------------------
    -- Storage arrays
    -- -------------------------------------------------------------------------
    type data_array   is array (0 to MAX_N-1)       of fp;
    type season_array is array (0 to MAX_M-1)       of fp;
    type savg_array   is array (0 to MAX_SEASONS-1) of fp;

    signal data_buf   : data_array   := (others => (others => '0'));
    signal S_arr      : season_array := (others => (others => '0'));
    signal season_avg : savg_array   := (others => (others => '0'));
    signal L_reg      : fp           := (others => '0');
    signal T_reg      : fp           := (others => '0');

    -- -------------------------------------------------------------------------
    -- FSM
    -- -------------------------------------------------------------------------
    type state_t is (
        S_IDLE,
        S_COLLECT,
        S_VALIDATE,
        S_WAIT_FC,
        S_INIT_L,
        S_INIT_T,
        S_INIT_SAVG,
        S_INIT_S,
        S_UPDATE,
        S_UPDATE_WAIT,
        S_FORECAST,
        S_ERROR
    );
    signal state : state_t := S_IDLE;

    -- -------------------------------------------------------------------------
    -- Counters
    -- -------------------------------------------------------------------------
    signal sample_cnt : integer range 0 to MAX_N       := 0;
    signal i_cnt      : integer range 0 to MAX_M       := 0;
    signal s_cnt      : integer range 0 to MAX_SEASONS := 0;
    signal t_cnt      : integer range 0 to MAX_N       := 0;
    signal k_cnt      : integer range 0 to MAX_HORIZON := 0;
    signal acc        : fp := (others => '0');

    -- UPDATE pipeline registers (2-stage to avoid combinational loop)
    signal upd_si    : integer range 0 to MAX_M-1 := 0;
    signal upd_fcast : fp := (others => '0');
    signal upd_S_cur : fp := (others => '0');
    signal upd_newL  : fp := (others => '0');
    signal upd_t     : integer range 0 to MAX_N-1 := 0;

    -- -------------------------------------------------------------------------
    -- Q16.16 multiply: a*b with rounding -> (a*b + 0.5) >> 16
    -- Product is 64-bit to avoid overflow
    -- -------------------------------------------------------------------------
    function fp_mul(a, b : fp) return fp is
        variable tmp : signed(63 downto 0);
    begin
        tmp := a * b + resize(to_signed(32768, 32), 64);
        return tmp(47 downto 16);
    end function;

    -- -------------------------------------------------------------------------
    -- Rounded Q16.16 reciprocal: 1/x  (integer x, result Q16.16)
    -- -------------------------------------------------------------------------
    function rcp(x : integer) return fp is
    begin
        return resize(to_signed((65536 + x/2) / x, 32), 32);
    end function;

    -- -------------------------------------------------------------------------
    -- Check if a divides b exactly
    -- -------------------------------------------------------------------------
    function divides(a, b : integer) return boolean is
    begin
        return (b mod a) = 0;
    end function;

begin

    process(clk)
        variable v_1ma       : fp;
        variable v_1mb       : fp;
        variable v_1mg       : fp;
        variable v_si        : integer range 0 to MAX_M-1;
        variable inv_M       : fp;
        variable inv_M2      : fp;
        variable inv_SEASONS : fp;
        variable idx_a       : integer range 0 to MAX_N-1;
        variable idx_b       : integer range 0 to MAX_N-1;
        variable k1_fp       : fp;
    begin
    if rising_edge(clk) then

        -- Default: pulse outputs are low unless explicitly driven
        forecast_valid <= '0';
        fitted_valid   <= '0';
        last_forecast  <= '0';
        error_out      <= '0';
        error_code     <= (others => '0');

        if rst = '1' then
            state      <= S_IDLE;
            sample_cnt <= 0;
            i_cnt      <= 0;
            s_cnt      <= 0;
            t_cnt      <= 0;
            k_cnt      <= 0;
            acc        <= (others => '0');
            L_reg      <= (others => '0');
            T_reg      <= (others => '0');

        else
            case state is

                -- =============================================================
                -- IDLE: wait for start pulse, latch M
                -- =============================================================
                when S_IDLE =>
                    if start = '1' then
                        M_act      <= m_in;
                        sample_cnt <= 0;
                        state      <= S_COLLECT;
                    end if;

                -- =============================================================
                -- COLLECT: count valid_in pulses into data_buf.
                --   - Collection ends when valid_in goes low (after >= 1 sample)
                --   - OR when buffer is full (MAX_N samples reached)
                -- =============================================================
                when S_COLLECT =>
                    if valid_in = '1' then
                        data_buf(sample_cnt) <= data_in;
                        if sample_cnt = MAX_N - 1 then
                            -- Buffer full
                            N_act      <= MAX_N;
                            sample_cnt <= 0;
                            acc        <= (others => '0');
                            i_cnt      <= 0;
                            state      <= S_VALIDATE;
                        else
                            sample_cnt <= sample_cnt + 1;
                        end if;
                    elsif sample_cnt > 0 then
                        -- valid_in deasserted after at least 1 sample received
                        N_act      <= sample_cnt;
                        sample_cnt <= 0;
                        acc        <= (others => '0');
                        i_cnt      <= 0;
                        state      <= S_VALIDATE;
                    end if;
                    -- sample_cnt=0 and valid_in=0: still waiting for first sample

                -- =============================================================
                -- VALIDATE: check N and M constraints
                -- =============================================================
                when S_VALIDATE =>
                    if M_act < 2 then
                        error_out  <= '1';
                        error_code <= "010";   -- M too small
                        state      <= S_ERROR;
                    elsif N_act < 2 * M_act then
                        error_out  <= '1';
                        error_code <= "001";   -- N < 2*M
                        state      <= S_ERROR;
                    elsif not divides(M_act, N_act) then
                        error_out  <= '1';
                        error_code <= "011";   -- M does not divide N
                        state      <= S_ERROR;
                    elsif (N_act / M_act) > MAX_SEASONS then
                        error_out  <= '1';
                        error_code <= "101";   -- too many seasons
                        state      <= S_ERROR;
                    else
                        SEASONS_act <= N_act / M_act;
                        state       <= S_WAIT_FC;
                    end if;

                -- =============================================================
                -- WAIT_FC: wait for forecast_start to latch horizon
                -- =============================================================
                when S_WAIT_FC =>
                    if forecast_start = '1' then
                        HORIZON_act <= horizon_in;
                        state       <= S_INIT_L;
                    end if;

                -- =============================================================
                -- INIT_L: L0 = mean(y[0..M-1])
                -- =============================================================
                when S_INIT_L =>
                    inv_M := rcp(M_act);
                    if i_cnt < M_act - 1 then
                        acc   <= acc + data_buf(i_cnt);
                        i_cnt <= i_cnt + 1;
                    else
                        L_reg <= fp_mul(acc + data_buf(M_act - 1), inv_M);
                        acc   <= (others => '0');
                        i_cnt <= 0;
                        state <= S_INIT_T;
                    end if;

                -- =============================================================
                -- INIT_T: T0 = sum(y[M+i] - y[i]) / M^2  for i in 0..M-1
                -- =============================================================
                when S_INIT_T =>
                    inv_M2 := rcp(M_act * M_act);
                    if i_cnt < M_act - 1 then
                        acc   <= acc + (data_buf(M_act + i_cnt) - data_buf(i_cnt));
                        i_cnt <= i_cnt + 1;
                    else
                        T_reg <= fp_mul(
                            acc + (data_buf(2*M_act - 1) - data_buf(M_act - 1)),
                            inv_M2
                        );
                        acc   <= (others => '0');
                        i_cnt <= 0;
                        s_cnt <= 0;
                        state <= S_INIT_SAVG;
                    end if;

                -- =============================================================
                -- INIT_SAVG: season_avg[s] = mean(y[s*M .. s*M+M-1])
                -- =============================================================
                when S_INIT_SAVG =>
                    inv_M := rcp(M_act);
                    idx_a := s_cnt * M_act + i_cnt;
                    idx_b := s_cnt * M_act + M_act - 1;
                    if i_cnt < M_act - 1 then
                        acc   <= acc + data_buf(idx_a);
                        i_cnt <= i_cnt + 1;
                    else
                        season_avg(s_cnt) <= fp_mul(acc + data_buf(idx_b), inv_M);
                        acc   <= (others => '0');
                        i_cnt <= 0;
                        if s_cnt = SEASONS_act - 1 then
                            s_cnt <= 0;
                            state <= S_INIT_S;
                        else
                            s_cnt <= s_cnt + 1;
                        end if;
                    end if;

                -- =============================================================
                -- INIT_S: S[i] = mean over seasons of (y[s*M+i] - savg[s])
                -- =============================================================
                when S_INIT_S =>
                    inv_SEASONS := rcp(SEASONS_act);
                    idx_a       := s_cnt * M_act + i_cnt;
                    idx_b       := (SEASONS_act - 1) * M_act + i_cnt;
                    if s_cnt < SEASONS_act - 1 then
                        acc   <= acc + (data_buf(idx_a) - season_avg(s_cnt));
                        s_cnt <= s_cnt + 1;
                    else
                        S_arr(i_cnt) <= fp_mul(
                            acc + (data_buf(idx_b) - season_avg(SEASONS_act - 1)),
                            inv_SEASONS
                        );
                        acc   <= (others => '0');
                        s_cnt <= 0;
                        if i_cnt = M_act - 1 then
                            i_cnt <= 0;
                            t_cnt <= 0;
                            state <= S_UPDATE;
                        else
                            i_cnt <= i_cnt + 1;
                        end if;
                    end if;

                -- =============================================================
                -- UPDATE stage 1: compute new values, latch into pipeline regs
                -- =============================================================
                when S_UPDATE =>
                    v_si      := t_cnt mod M_act;
                    v_1ma     := ONE - alpha;
                    upd_fcast <= L_reg + T_reg + S_arr(v_si);
                    upd_newL  <= fp_mul(alpha, data_buf(t_cnt) - S_arr(v_si))
                               + fp_mul(v_1ma, L_reg + T_reg);
                    upd_S_cur <= S_arr(v_si);
                    upd_si    <= v_si;
                    upd_t     <= t_cnt;
                    state     <= S_UPDATE_WAIT;

                -- =============================================================
                -- UPDATE stage 2: commit L, T, S; output fitted value
                -- =============================================================
                when S_UPDATE_WAIT =>
                    v_1mb := ONE - beta;
                    v_1mg := ONE - gamma;

                    -- Output fitted (one-step-ahead prediction before update)
                    fitted_out   <= upd_fcast;
                    fitted_valid <= '1';

                    -- Commit updated state
                    L_reg         <= upd_newL;
                    T_reg         <= fp_mul(beta,  upd_newL - L_reg)
                                   + fp_mul(v_1mb, T_reg);
                    S_arr(upd_si) <= fp_mul(gamma, data_buf(upd_t) - upd_newL)
                                   + fp_mul(v_1mg, upd_S_cur);

                    if t_cnt = N_act - 1 then
                        k_cnt <= 0;
                        state <= S_FORECAST;
                    else
                        t_cnt <= t_cnt + 1;
                        state <= S_UPDATE;
                    end if;

                -- =============================================================
                -- FORECAST: output HORIZON_act forecasts
                --   F(k) = L + k*T + S[(N + k - 1) mod M]   k = 1..HORIZON
                -- last_forecast pulses on the final output word
                -- =============================================================
                when S_FORECAST =>
                    k1_fp := to_signed((k_cnt + 1) * 65536, 32);

                    forecast_out   <= L_reg
                                    + fp_mul(k1_fp, T_reg)
                                    + S_arr((N_act + k_cnt) mod M_act);
                    forecast_valid <= '1';

                    if k_cnt = HORIZON_act - 1 then
                        last_forecast <= '1';
                        state         <= S_IDLE;
                    else
                        k_cnt <= k_cnt + 1;
                    end if;

                -- =============================================================
                -- ERROR: hold error_out for one cycle then return to idle
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
