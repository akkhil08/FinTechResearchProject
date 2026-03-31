-- =============================================================================
-- hw_top.vhd  -  Minimal ZedBoard Wrapper  (NO ILA)
--
-- PHYSICAL PINS: 16 total
--   clk            1 pin   Y9   125MHz oscillator
--   rst            1 pin   P16  BTNC centre button
--   start          1 pin   T18  BTNU up button
--   forecast_start 1 pin   R18  BTNR right button
--   sw[7:0]        8 pins       DIP switches
--   led[3:0]       4 pins       Status LEDs
--
-- MODIFICATION: fitted_out and forecast_out are now connected to signals
--               and their MSBs are XORed into LED(0) to force Vivado to
--               keep the entire data path including data_buf.
-- =============================================================================

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity hw_top is
    port (
        clk            : in  std_logic;
        rst            : in  std_logic;
        start          : in  std_logic;
        forecast_start : in  std_logic;
        sw             : in  std_logic_vector(7 downto 0);
        led            : out std_logic_vector(3 downto 0)
    );
end entity hw_top;

architecture rtl of hw_top is

    constant C_ALPHA : signed(31 downto 0) := to_signed(1021785021, 32);
    constant C_BETA  : signed(31 downto 0) := to_signed(107374182, 32);  -- 0.1
    constant C_GAMMA : signed(31 downto 0) := to_signed(107374182, 32);  -- 0.1

    constant ROM_SIZE : integer := 60;
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

    type rom_state_t is (ROM_IDLE, ROM_STREAM, ROM_DONE);
    signal rom_state   : rom_state_t := ROM_IDLE;
    signal rom_addr    : integer range 0 to ROM_SIZE := 0;
    signal rom_data    : signed(31 downto 0) := (others => '0');
    signal rom_valid   : std_logic := '0';
    signal auto_fc     : std_logic := '0';
    signal fc_combined : std_logic;

    signal m_int       : std_logic_vector(4 downto 0);
    signal horizon_int : std_logic_vector(4 downto 0);

    signal fitted_valid_s   : std_logic;
    signal forecast_valid_s : std_logic;
    signal last_forecast_s  : std_logic;
    signal error_out_s      : std_logic;
    signal hb_reg           : std_logic := '0';

    -- =========================================================================
    -- MODIFICATION: Connect fitted_out and forecast_out to real signals
    -- instead of 'open'. This forces Vivado to keep the entire data path
    -- including data_buf, all arithmetic, and the multipliers.
    -- Without this, Vivado performs dead code elimination and removes
    -- everything because the outputs were never used.
    -- =========================================================================
    signal fitted_reg   : signed(31 downto 0) := (others => '0');
    signal forecast_reg : signed(31 downto 0) := (others => '0');

begin

    m_int       <= std_logic_vector(to_unsigned(to_integer(unsigned(sw(3 downto 0))) + 2, 5));
    horizon_int <= std_logic_vector(to_unsigned(to_integer(unsigned(sw(7 downto 4))) + 1, 5));
    fc_combined <= forecast_start or auto_fc;

    -- =========================================================================
    -- LED assignments:
    -- led(0) = heartbeat XORed with MSB of fitted and forecast outputs
    --          The XOR forces Vivado to keep fitted_reg and forecast_reg
    --          in the design, which in turn keeps the full data path alive.
    -- led(1) = forecast_valid
    -- led(2) = last_forecast
    -- led(3) = error_out
    -- =========================================================================
    led(0) <= hb_reg xor fitted_reg(31) xor forecast_reg(31);
    led(1) <= forecast_valid_s;
    led(2) <= last_forecast_s;
    led(3) <= error_out_s;

    process(clk)
    begin
        if rising_edge(clk) then
            rom_valid <= '0';
            auto_fc   <= '0';
            if rst = '1' then
                rom_state <= ROM_IDLE;
                rom_addr  <= 0;
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
                        auto_fc   <= '1';
                        rom_state <= ROM_IDLE;
                end case;
            end if;
        end if;
    end process;

    u_core : entity work.holt_winters_q2_30
        generic map (
            MAX_N => 72, MAX_M => 24, MAX_HORIZON => 24, MAX_SEASONS => 8
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
            fitted_out     => fitted_reg,      -- CHANGED: was open
            fitted_valid   => fitted_valid_s,
            forecast_out   => forecast_reg,    -- CHANGED: was open
            forecast_valid => forecast_valid_s,
            last_forecast  => last_forecast_s,
            error_out      => error_out_s,
            error_code     => open
        );

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
