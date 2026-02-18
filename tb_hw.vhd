library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
-- use std.env.all; -- optional, replaced with assert for Vivado-safe finish

entity tb_holt_winters is
end tb_holt_winters;

architecture sim of tb_holt_winters is

    signal clk, rst, start : std_logic := '0';
    signal valid_in : std_logic := '0';
    signal data_in  : signed(31 downto 0);

    signal alpha, beta, gamma : signed(31 downto 0);
    signal forecast_out   : signed(31 downto 0);
    signal forecast_valid : std_logic;

    constant clk_period : time := 10 ns;

    type price_array is array (0 to 59) of integer;
    constant PSEI_2017_2021 : price_array := (
        7229, 7212, 7311, 7661, 7837, 8018, 7977, 7954, 7853, 7970, 8029, 8053,
        8022, 7977, 7953, 8029, 8054, 8022, 7977, 7953, 8029, 8054, 8022, 7977,
        7953, 8029, 8054, 8022, 7977, 7953, 8029, 8054, 8022, 7977, 7953, 8029,
        8054, 8022, 7201, 7055, 6953, 6855, 6270, 6902, 6628, 6371, 6443, 6795,
        6613, 7139, 6791, 6324, 5864, 5884, 5928, 6202, 5812, 5212, 5732, 7122
    );

    type forecast_array is array (0 to 11) of integer;
    signal forecasts : forecast_array := (others => 0);

begin

    -- Clock generation
    clk <= not clk after clk_period / 2;

    -- DUT instantiation
    dut : entity work.holt_winters
        port map (
            clk => clk, rst => rst, start => start,
            alpha => alpha, beta => beta, gamma => gamma,
            data_in => data_in, valid_in => valid_in,
            forecast_out => forecast_out, forecast_valid => forecast_valid
        );

    -- Stimulus process
    stim_proc : process
        variable price_q16 : signed(31 downto 0);
        variable i : integer;
    begin
        -- Reset
        rst <= '1'; wait for 20 ns; rst <= '0';

        -- Q16.16 coefficients (alpha=0.95, beta=0, gamma=0)
        alpha <= to_signed(integer(0.951593*65536.0), 32);
        beta  <= to_signed(0, 32);
        gamma <= to_signed(0, 32);

        -- Start signal aligned to clock
        start <= '1';
        wait until rising_edge(clk);
        start <= '0';

        -- Feed input data
        i := 0;
        while i <= 59 loop
            price_q16 := to_signed(PSEI_2017_2021(i), 32) sll 16;  -- safer Q16.16
            data_in <= price_q16;
            valid_in <= '1';
            wait until rising_edge(clk);
            i := i + 1;
        end loop;

        valid_in <= '0';
        wait until rising_edge(clk);

        -- Wait for forecasts to complete
        wait for 300 ns;
        -- End simulation safely
        assert false report "Simulation Finished" severity failure;
    end process;

    -- Monitor process
    monitor_proc : process
        variable count : integer := 0;
    begin
        loop
            wait until rising_edge(clk);

            if forecast_valid = '1' then
                report "Forecast Month " & integer'image(count + 1) & 
                       ": " & integer'image(to_integer(forecast_out srl 16))
                       severity note;

                count := count + 1;
                exit when count = 12;
            end if;
        end loop;

        report "FINAL 12 FORECASTS COMPLETE" severity note;
        wait for 100 ns;
        -- End simulation safely
        assert false report "Simulation Finished" severity failure;
    end process;

end sim;
