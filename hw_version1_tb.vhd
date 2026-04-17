library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity holt_winters_tb is
end holt_winters_tb;

architecture sim of holt_winters_tb is

    signal clk            : std_logic := '0';
    signal rst            : std_logic := '0';
    signal start          : std_logic := '0';
    signal alpha          : signed(31 downto 0) := (others => '0');
    signal beta           : signed(31 downto 0) := (others => '0');
    signal gamma          : signed(31 downto 0) := (others => '0');
    signal data_in        : signed(31 downto 0) := (others => '0');
    signal valid_in       : std_logic := '0';
    signal fitted_out     : signed(31 downto 0);
    signal fitted_valid   : std_logic;
    signal forecast_out   : signed(31 downto 0);
    signal forecast_valid : std_logic;
    signal done           : std_logic;

    constant clk_period : time := 10 ns;

    type price_array is array (0 to 59) of signed(31 downto 0);
    constant PSEI_DATA : price_array := (
        to_signed(473802998,32), to_signed(472651530,32), to_signed(479180882,32),
        to_signed(502071951,32), to_signed(513613496,32), to_signed(514009334,32),
        to_signed(525470925,32), to_signed(521572844,32), to_signed(535522836,32),
        to_signed(548225679,32), to_signed(540936110,32), to_signed(560884613,32),
        to_signed(574358159,32), to_signed(555436605,32), to_signed(522966139,32),
        to_signed(512442368,32), to_signed(491334533,32), to_signed(471445012,32),
        to_signed(502792192,32), to_signed(514831811,32), to_signed(476893676,32),
        to_signed(467946045,32), to_signed(482859418,32), to_signed(489293087,32),
        to_signed(524778209,32), to_signed(504986993,32), to_signed(519106068,32),
        to_signed(521189458,32), to_signed(522323231,32), to_signed(524268995,32),
        to_signed(527289549,32), to_signed(522954998,32), to_signed(509809132,32),
        to_signed(522788536,32), to_signed(507180483,32), to_signed(512180879,32),
        to_signed(471910973,32), to_signed(444852470,32), to_signed(348732129,32),
        to_signed(373601731,32), to_signed(382654218,32), to_signed(406829138,32),
        to_signed(388526899,32), to_signed(385625620,32), to_signed(384318177,32),
        to_signed(414449664,32), to_signed(445085123,32), to_signed(467908035,32),
        to_signed(433364664,32), to_signed(445307945,32), to_signed(422254346,32),
        to_signed(417521336,32), to_signed(434404721,32), to_signed(452323574,32),
        to_signed(410925793,32), to_signed(449278116,32), to_signed(455663944,32),
        to_signed(462336819,32), to_signed(471916872,32), to_signed(466788680,32)
    );

    -- -------------------------------------------------------------------------
    -- Convert Q16.16 signed value to "XXXXX.XX" decimal string
    -- Integer part  = raw / 65536
    -- Fractional 2dp= ((raw mod 65536) * 100) / 65536
    -- -------------------------------------------------------------------------
    function q16_to_str(q : signed(31 downto 0)) return string is
        variable raw  : integer;
        variable ipart: integer;
        variable fpart: integer;
        variable neg  : boolean;
    begin
        raw  := to_integer(q);
        neg  := raw < 0;
        if neg then raw := -raw; end if;
        ipart := raw / 65536;
        fpart := ((raw mod 65536) * 100) / 65536;
        if neg then
            if    fpart = 0  then return "-" & integer'image(ipart) & ".00";
            elsif fpart < 10 then return "-" & integer'image(ipart) & ".0" & integer'image(fpart);
            else               return "-" & integer'image(ipart) & "." & integer'image(fpart);
            end if;
        else
            if    fpart = 0  then return integer'image(ipart) & ".00";
            elsif fpart < 10 then return integer'image(ipart) & ".0" & integer'image(fpart);
            else               return integer'image(ipart) & "." & integer'image(fpart);
            end if;
        end if;
    end function;

begin

    clk <= not clk after clk_period / 2;

    dut : entity work.holt_winters
        port map (
            clk=>clk, rst=>rst, start=>start,
            alpha=>alpha, beta=>beta, gamma=>gamma,
            data_in=>data_in, valid_in=>valid_in,
            fitted_out=>fitted_out, fitted_valid=>fitted_valid,
            forecast_out=>forecast_out, forecast_valid=>forecast_valid,
            done=>done
        );

    -- =========================================================================
    -- STIM
    -- =========================================================================
    stim_proc : process
    begin
        report "SIM START" severity note;

        rst <= '1'; wait for 20 ns; rst <= '0'; wait for 20 ns;
        report "RESET DONE" severity note;

        alpha <= to_signed(62364, 32);
        beta  <= to_signed(0, 32);
        gamma <= to_signed(0, 32);

        wait until rising_edge(clk);
        start <= '1';
        wait until rising_edge(clk);
        start <= '0';
        report "START PULSED" severity note;

        for i in 0 to 59 loop
            data_in  <= PSEI_DATA(i);
            valid_in <= '1';
            wait until rising_edge(clk);
        end loop;
        valid_in <= '0';
        report "ALL 60 SAMPLES SENT" severity note;

        wait for 10000 ns;

        report "TIMEOUT - done never came" severity note;
        assert false report "SIM END" severity failure;
    end process;

    -- =========================================================================
    -- DONE MONITOR
    -- =========================================================================
    done_mon : process
    begin
        wait until done = '1';
        report "DONE SIGNAL RECEIVED" severity note;
        wait;
    end process;

    -- =========================================================================
    -- FORECAST MONITOR
    -- Uses q16_to_str to display decimal places e.g. "7062.11" not "7062"
    -- =========================================================================
    monitor_proc : process
        variable count : integer := 0;
    begin
        loop
            wait until rising_edge(clk);
            if forecast_valid = '1' then
                report "Forecast Month " & integer'image(count + 1) &
                       " : PSEi = " & q16_to_str(forecast_out)
                    severity note;
                count := count + 1;
                exit when count = 12;
            end if;
        end loop;
        report "ALL 12 FORECASTS DONE" severity note;
        assert false report "Monitor done" severity failure;
    end process;

end sim;
