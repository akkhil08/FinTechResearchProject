-- =============================================================================
-- hw_top.vhd  -  Minimal ZedBoard Wrapper (v10)
--
-- Updated for holt_winters_q2_30_opt_v10:
--   MAX_SEASONS now carries the actual N/M value computed by the PS.
--   For this config: 60 samples / 12 per season = 5 seasons.
--   No separate INIT_SEASONS generic needed -- MAX_SEASONS does the job.
-- =============================================================================

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity hw_top is
    port (
        clk   : in  std_logic;
        rst   : in  std_logic;
        start : in  std_logic;
        sw    : in  std_logic_vector(7 downto 0);
        led   : out std_logic_vector(3 downto 0)
    );
end entity hw_top;

architecture rtl of hw_top is

    constant C_ALPHA : signed(31 downto 0) := to_signed(1021785021, 32);
    constant C_BETA  : signed(31 downto 0) := to_signed(107374182,  32);
    constant C_GAMMA : signed(31 downto 0) := to_signed(107374182,  32);

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

    type rom_state_t is (ROM_IDLE, ROM_PREFETCH, ROM_STREAM, ROM_DONE);
    signal rom_state : rom_state_t := ROM_IDLE;
    signal rom_addr  : integer range 0 to ROM_SIZE := 0;
    signal rom_data  : signed(31 downto 0) := (others => '0');
    signal rom_valid : std_logic := '0';

    signal fitted_valid_s   : std_logic;
    signal forecast_valid_s : std_logic;
    signal last_forecast_s  : std_logic;
    signal error_out_s      : std_logic;
    signal hb_reg           : std_logic := '0';

    signal fitted_reg   : signed(31 downto 0) := (others => '0');
    signal forecast_reg : signed(31 downto 0) := (others => '0');

begin

    led(0) <= hb_reg xor fitted_reg(31) xor forecast_reg(31);
    led(1) <= forecast_valid_s;
    led(2) <= last_forecast_s;
    led(3) <= error_out_s;

    -- =========================================================================
    -- ROM streamer (unchanged from v9)
    -- =========================================================================
    process(clk)
    begin
        if rising_edge(clk) then
            rom_valid <= '0';
            if rst = '1' then
                rom_state <= ROM_IDLE;
                rom_addr  <= 0;
                rom_data  <= (others => '0');
            else
                case rom_state is
                    when ROM_IDLE =>
                        if start = '1' then
                            rom_addr  <= 0;
                            rom_state <= ROM_PREFETCH;
                        end if;
                    when ROM_PREFETCH =>
                        rom_data  <= PSEI_ROM(rom_addr);
                        rom_addr  <= rom_addr + 1;
                        rom_state <= ROM_STREAM;
                    when ROM_STREAM =>
                        rom_valid <= '1';
                        if rom_addr = ROM_SIZE then
                            rom_state <= ROM_DONE;
                        else
                            rom_data <= PSEI_ROM(rom_addr);
                            rom_addr <= rom_addr + 1;
                        end if;
                    when ROM_DONE =>
                        rom_addr  <= 0;
                        rom_state <= ROM_IDLE;
                end case;
            end if;
        end if;
    end process;

    -- =========================================================================
    -- Core instantiation.
    --
    -- MAX_SEASONS is set to N/M = 60/12 = 5, computed by the PS.
    -- The core uses MAX_SEASONS directly as SEASONS_act -- no PL-side divide.
    --
    -- Rule: always keep MAX_SEASONS = MAX_N / MAX_M exactly.
    -- =========================================================================
    u_core : entity work.holt_winters_q2_30
        generic map (
            MAX_N       => 60,
            MAX_M       => 12,
            MAX_HORIZON => 12,
            MAX_SEASONS => 5    -- PS computes: 60 / 12 = 5
        )
        port map (
            clk            => clk,
            rst            => rst,
            alpha          => C_ALPHA,
            beta           => C_BETA,
            gamma          => C_GAMMA,
            data_in        => rom_data,
            valid_in       => rom_valid,
            fitted_out     => fitted_reg,
            fitted_valid   => fitted_valid_s,
            forecast_out   => forecast_reg,
            forecast_valid => forecast_valid_s,
            last_forecast  => last_forecast_s,
            error_out      => error_out_s,
            error_code     => open
        );

    -- =========================================================================
    -- Heartbeat
    -- =========================================================================
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