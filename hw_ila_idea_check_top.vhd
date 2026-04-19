-- =============================================================================
-- hw_top.vhd  -  ZedBoard Hardware-Test Wrapper + Internal ROM + System ILA
--
-- NO EXTERNAL HARDWARE NEEDED. No jumper wires. No UART dongle.
-- The 60 PSEi samples are embedded as a ROM inside this wrapper.
--
-- COMPLETE PROCEDURE (only 2 buttons):
--   1. Set DIP switches FIRST:
--        SW[3:0] = 1010  ->  M = 12  (monthly seasons)
--        SW[7:4] = 1011  ->  H = 12  (12-month forecast)
--   2. In Vivado ILA: set trigger on probe8==1, click "Run Trigger"
--   3. Press BTNC (Center) = rst
--   4. Press BTNU (Up)     = start
--      -> ROM streams all 60 samples automatically (one per clock)
--      -> FSM collects, validates, initialises
--      -> forecast_start fires automatically after last sample
--      -> 12 forecasts output
--      -> ILA captures fitted values + forecasts
--
-- BTNR (forecast_start) kept as manual override but not required.
-- =============================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity holt_winters_top is
    port (
        clk             : in  std_logic;
        rst             : in  std_logic;
        start           : in  std_logic;
        forecast_start  : in  std_logic;
        sw              : in  std_logic_vector(7 downto 0);
        fitted_valid    : out std_logic;
        forecast_valid  : out std_logic;
        last_forecast   : out std_logic;
        error_out       : out std_logic;
        error_code      : out std_logic_vector(2 downto 0);
        heartbeat       : out std_logic;
        forecast_byte   : out std_logic_vector(7 downto 0)
    );
end entity holt_winters_top;

architecture rtl of holt_winters_top is

    -- -------------------------------------------------------------------------
    -- Smoothing constants  alpha=0.95, beta=0, gamma=0  (Q16.16)
    -- -------------------------------------------------------------------------
    constant C_ALPHA : signed(31 downto 0) := to_signed(62364, 32);
    constant C_BETA  : signed(31 downto 0) := to_signed(0,     32);
    constant C_GAMMA : signed(31 downto 0) := to_signed(0,     32);

    constant ROM_SIZE : integer := 60;

    -- -------------------------------------------------------------------------
    -- 60-sample PSEi ROM  (Q16.16 fixed-point, real = entry/65536)
    -- -------------------------------------------------------------------------
    type rom_t is array (0 to ROM_SIZE-1) of signed(31 downto 0);
    constant PSEI_ROM : rom_t := (
        to_signed(473802998, 32), to_signed(472651530, 32), to_signed(479180882, 32),
        to_signed(502071951, 32), to_signed(513613496, 32), to_signed(514009334, 32),
        to_signed(525470925, 32), to_signed(521572844, 32), to_signed(535522836, 32),
        to_signed(548225679, 32), to_signed(540936110, 32), to_signed(560884613, 32),
        to_signed(574358159, 32), to_signed(555436605, 32), to_signed(522966139, 32),
        to_signed(512442368, 32), to_signed(491334533, 32), to_signed(471445012, 32),
        to_signed(502792192, 32), to_signed(514831811, 32), to_signed(476893676, 32),
        to_signed(467946045, 32), to_signed(482859418, 32), to_signed(489293087, 32),
        to_signed(524778209, 32), to_signed(504986993, 32), to_signed(519106068, 32),
        to_signed(521189458, 32), to_signed(522323231, 32), to_signed(524268995, 32),
        to_signed(527289549, 32), to_signed(522954998, 32), to_signed(509809132, 32),
        to_signed(522788536, 32), to_signed(507180483, 32), to_signed(512180879, 32),
        to_signed(471910973, 32), to_signed(444852470, 32), to_signed(348732129, 32),
        to_signed(373601731, 32), to_signed(382654218, 32), to_signed(406829138, 32),
        to_signed(388526899, 32), to_signed(385625620, 32), to_signed(384318177, 32),
        to_signed(414449664, 32), to_signed(445085123, 32), to_signed(467908035, 32),
        to_signed(433364664, 32), to_signed(445307945, 32), to_signed(422254346, 32),
        to_signed(417521336, 32), to_signed(434404721, 32), to_signed(452323574, 32),
        to_signed(410925793, 32), to_signed(449278116, 32), to_signed(455663944, 32),
        to_signed(462336819, 32), to_signed(471916872, 32), to_signed(466788680, 32)
    );

    -- -------------------------------------------------------------------------
    -- ROM streamer signals
    -- -------------------------------------------------------------------------
    type rom_state_t is (ROM_IDLE, ROM_STREAM, ROM_DONE);
    signal rom_state : rom_state_t := ROM_IDLE;
    signal rom_addr  : integer range 0 to ROM_SIZE := 0;
    signal rom_data  : signed(31 downto 0) := (others => '0');
    signal rom_valid : std_logic := '0';
    signal auto_fc   : std_logic := '0';
    signal fc_combined : std_logic;

    -- -------------------------------------------------------------------------
    -- DUT signals
    -- -------------------------------------------------------------------------
    signal m_int            : integer range 2 to 24;
    signal horizon_int      : integer range 1 to 24;
    signal fitted_out_s     : signed(31 downto 0);
    signal fitted_valid_s   : std_logic;
    signal forecast_out_s   : signed(31 downto 0);
    signal forecast_valid_s : std_logic;
    signal last_forecast_s  : std_logic;
    signal error_out_s      : std_logic;
    signal error_code_s     : std_logic_vector(2 downto 0);
    signal hb_reg           : std_logic := '0';

    -- ILA probe wrappers
    signal pr1 : std_logic_vector(0 downto 0);
    signal pr3 : std_logic_vector(0 downto 0);
    signal pr4 : std_logic_vector(0 downto 0);
    signal pr5 : std_logic_vector(0 downto 0);
    signal pr7 : std_logic_vector(0 downto 0);
    signal pr8 : std_logic_vector(0 downto 0);
    signal pr9 : std_logic_vector(0 downto 0);

    component ila_holt_winters
        port (
            clk    : in std_logic;
            probe0 : in std_logic_vector(31 downto 0);
            probe1 : in std_logic_vector(0  downto 0);
            probe2 : in std_logic_vector(31 downto 0);
            probe3 : in std_logic_vector(0  downto 0);
            probe4 : in std_logic_vector(0  downto 0);
            probe5 : in std_logic_vector(0  downto 0);
            probe6 : in std_logic_vector(2  downto 0);
            probe7 : in std_logic_vector(0  downto 0);
            probe8 : in std_logic_vector(0  downto 0);
            probe9 : in std_logic_vector(0  downto 0)
        );
    end component;

begin

    -- Output assignments
    fitted_valid   <= fitted_valid_s;
    forecast_valid <= forecast_valid_s;
    last_forecast  <= last_forecast_s;
    error_out      <= error_out_s;
    error_code     <= error_code_s;
    heartbeat      <= hb_reg;
    forecast_byte  <= std_logic_vector(forecast_out_s(23 downto 16));

    -- forecast_start = manual button OR auto after ROM finishes
    fc_combined <= forecast_start or auto_fc;

    -- SW decode
    m_int       <= to_integer(unsigned(sw(3 downto 0))) + 2;
    horizon_int <= to_integer(unsigned(sw(7 downto 4))) + 1;

    -- ILA probes
    pr1(0) <= fitted_valid_s;
    pr3(0) <= forecast_valid_s;
    pr4(0) <= last_forecast_s;
    pr5(0) <= error_out_s;
    pr7(0) <= rom_valid;
    pr8(0) <= start;
    pr9(0) <= hb_reg;

    -- -------------------------------------------------------------------------
    -- ROM streamer: triggered by start pulse
    -- Streams one sample per clock, then auto-fires forecast_start
    -- -------------------------------------------------------------------------
    process(clk)
    begin
        if rising_edge(clk) then
            rom_valid <= '0';
            auto_fc   <= '0';

            if rst = '1' then
                rom_state <= ROM_IDLE;
                rom_addr  <= 0;
                rom_data  <= (others => '0');

            else
                case rom_state is

                    when ROM_IDLE =>
                        if start = '1' then
                            rom_addr  <= 0;
                            rom_state <= ROM_STREAM;
                        end if;

                    when ROM_STREAM =>
                        rom_data  <= PSEI_ROM(rom_addr);
                        rom_valid <= '1';
                        if rom_addr = ROM_SIZE - 1 then
                            rom_state <= ROM_DONE;
                        else
                            rom_addr <= rom_addr + 1;
                        end if;

                    when ROM_DONE =>
                        auto_fc   <= '1';   -- pulse forecast_start for 1 cycle
                        rom_state <= ROM_IDLE;

                end case;
            end if;
        end if;
    end process;

    -- -------------------------------------------------------------------------
    -- Holt-Winters DUT
    -- -------------------------------------------------------------------------
    u_dut : entity work.holt_winters
        generic map (
            MAX_N       => 72,
            MAX_M       => 24,
            MAX_HORIZON => 24,
            MAX_SEASONS => 8
        )
        port map (
            clk            => clk,
            rst            => rst,
            m_in           => m_int,
            forecast_start => fc_combined,
            horizon_in     => horizon_int,
            start          => start,
            alpha          => C_ALPHA,
            beta           => C_BETA,
            gamma          => C_GAMMA,
            data_in        => rom_data,
            valid_in       => rom_valid,
            fitted_out     => fitted_out_s,
            fitted_valid   => fitted_valid_s,
            forecast_out   => forecast_out_s,
            forecast_valid => forecast_valid_s,
            last_forecast  => last_forecast_s,
            error_out      => error_out_s,
            error_code     => error_code_s
        );

    -- -------------------------------------------------------------------------
    -- ILA
    -- -------------------------------------------------------------------------
    u_ila : ila_holt_winters
        port map (
            clk    => clk,
            probe0 => std_logic_vector(fitted_out_s),
            probe1 => pr1,
            probe2 => std_logic_vector(forecast_out_s),
            probe3 => pr3,
            probe4 => pr4,
            probe5 => pr5,
            probe6 => error_code_s,
            probe7 => pr7,
            probe8 => pr8,
            probe9 => pr9
        );

    -- -------------------------------------------------------------------------
    -- Heartbeat: toggles on every fitted_valid pulse
    -- -------------------------------------------------------------------------
    process(clk)
    begin
        if rising_edge(clk) then
            if rst = '1' then
                hb_reg <= '0';
            elsif fitted_valid_s = '1' then
                hb_reg <= not hb_reg;
            end if;
        end if;
    end process;

end architecture rtl;
