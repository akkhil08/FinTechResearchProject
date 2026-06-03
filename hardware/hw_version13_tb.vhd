-- =============================================================================
-- hw_tb_hardcore_v17.vhd
--
-- Testbench for holt_winters_q2_30_opt_v17.
-- CLK_PERIOD = 10 ns (100 MHz). Single test: beta=0, gamma=0.
--
-- CHANGES FROM v11 TESTBENCH:
--   1. CLK_PERIOD: 14 ns -> 10 ns  (100 MHz target)
--   2. CLK_NS:     14    -> 10     (used in latency calculations)
--   3. Version label updated in report strings
--   4. Expected cycle counts will be higher due to deeper pipeline in v17
--      (v17 has 18 FSM states vs 13 in v11, ~649 cycles vs ~355 cycles)
--      All cycle measurements are DYNAMIC -- no hardcoded expected counts.
--   5. Streaming protocol UNCHANGED -- valid_in only, same 3-step protocol.
--   6. Entity ports UNCHANGED -- same port map as v11.
--   7. Expected forecast VALUES unchanged -- same algorithm, same dataset.
--   8. FC_TOLERANCE unchanged -- algorithm correctness independent of timing.
--
-- STREAMING PROTOCOL (unchanged from v11):
--   STEP 1: data_in=ROM(0), valid_in='0'  -- prefetch, data stable
--   STEP 2: data_in=ROM(0), valid_in='1'  -- S_IDLE fires, goes to S_COLLECT
--   STEP 3: data_in=ROM(0), valid_in='1'  -- S_COLLECT writes ROM(0) to buf(0)
--   STEP 4: data_in=ROM(1..59), valid_in='1' -- stream remaining samples
--
--   Total valid_in='1' pulses: 61
--   Total samples written: 60 (data_buf[0..59])
--
-- PIPELINE DEPTH NOTE (v17 vs v11):
--   v11: 13 states, ~355 cycles total, UPDATE = 4 cycles/sample
--   v17: 18 states, ~649 cycles total, UPDATE = 8 cycles/sample
--   Extra states: S_INIT_LT_COMMIT, S_UPDATE_MULT2, S_UPDATE_WAIT2,
--                 S_UPDATE_FINALIZE2, S_FORECAST_PRE (+5 vs v11)
--   All extra states are pipeline registers to close 100 MHz timing.
--   Algorithmic output (forecast values) is identical to v11.
-- =============================================================================

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity tb_holt_winters_q2_30_timing is
end entity;

architecture sim of tb_holt_winters_q2_30_timing is

    -- =========================================================================
    -- CHANGED: 14 ns -> 10 ns for 100 MHz
    -- =========================================================================
    constant CLK_PERIOD : time    := 10 ns;
    constant TSETUP     : time    :=  1 ns;
    constant CLK_NS     : integer := 10;  -- CHANGED: 14 -> 10

    constant C_ALPHA : integer := 1021785021;

    signal clk : std_logic := '0';
    signal rst : std_logic := '1';

    signal alpha : signed(31 downto 0) := to_signed(C_ALPHA, 32);
    signal beta  : signed(31 downto 0) := to_signed(0, 32);
    signal gamma : signed(31 downto 0) := to_signed(0, 32);

    signal data_in        : signed(31 downto 0) := (others => '0');
    signal valid_in       : std_logic := '0';

    signal fitted_out     : signed(31 downto 0);
    signal fitted_valid   : std_logic;
    signal forecast_out   : signed(31 downto 0);
    signal forecast_valid : std_logic;
    signal last_forecast  : std_logic;
    signal error_out      : std_logic;
    signal error_code     : std_logic_vector(2 downto 0);

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

    -- =========================================================================
    -- Expected forecast values (beta=0, gamma=0, correct 60-sample dataset)
    -- UNCHANGED from v11 -- same algorithm produces same outputs
    -- =========================================================================
    type fc_expected_t is array (0 to 11) of integer;
    constant EXPECTED_FC : fc_expected_t := (
        462798048, 450932087, 423863058, 429910177,
        432567222, 436579428, 432973611, 439916351,
        432635091, 442472491, 448048411, 456993486
    );
    constant FC_TOLERANCE : integer := 6553600; -- ~100 in Q16.16

    signal cycle_count     : integer := 0;
    signal abs_cycle_count : integer := 0;

    shared variable sv_start_snap     : integer := 0;
    shared variable sv_start_snap_set : boolean := false;

    signal first_fitted_cycle    : integer := 0;
    signal second_fitted_cycle   : integer := 0;
    signal last_fitted_cycle     : integer := 0;
    signal first_forecast_cycle  : integer := 0;
    signal second_forecast_cycle : integer := 0;
    signal last_forecast_cycle   : integer := 0;

    signal fitted_count   : integer := 0;
    signal forecast_count : integer := 0;

    signal first_fitted_seen    : boolean := false;
    signal second_fitted_seen   : boolean := false;
    signal first_forecast_seen  : boolean := false;
    signal second_forecast_seen : boolean := false;

    signal done : boolean := false;

    constant MAX_HORIZON_TB : integer := 24;
    type int_array is array (0 to MAX_HORIZON_TB-1) of integer;

    shared variable sv_fc_log_ipart : int_array := (others => 0);
    shared variable sv_fc_log_fpart : int_array := (others => 0);
    shared variable sv_fc_log_cycle : int_array := (others => 0);
    shared variable sv_fc_log_raw   : int_array := (others => 0);

    procedure decode_q1616(
        raw   : in  integer;
        ipart : out integer;
        fpart : out integer
    ) is
        variable ip : integer;
        variable fr : integer;
    begin
        ip := raw / 65536;
        fr := raw - ip * 65536;
        if fr < 0 then
            ip := ip - 1;
            fr := fr + 65536;
        end if;
        ipart := ip;
        fpart := fr * 10000 / 65536;
    end procedure;

begin

    clk <= not clk after CLK_PERIOD / 2;

    -- =========================================================================
    -- DUT instantiation (port map identical to v11 testbench)
    -- =========================================================================
    u_dut : entity work.holt_winters_q2_30
        generic map (
            MAX_N       => 60,
            MAX_M       => 12,
            MAX_HORIZON => 12,
            MAX_SEASONS => 5
        )
        port map (
            clk            => clk,
            rst            => rst,
            alpha          => alpha,
            beta           => beta,
            gamma          => gamma,
            data_in        => data_in,
            valid_in       => valid_in,
            fitted_out     => fitted_out,
            fitted_valid   => fitted_valid,
            forecast_out   => forecast_out,
            forecast_valid => forecast_valid,
            last_forecast  => last_forecast,
            error_out      => error_out,
            error_code     => error_code
        );

    -- =========================================================================
    -- Cycle counter
    -- =========================================================================
    process(clk)
    begin
        if rising_edge(clk) then
            abs_cycle_count <= abs_cycle_count + 1;
            cycle_count     <= cycle_count + 1;
        end if;
    end process;

    -- =========================================================================
    -- Output monitor (logic identical to v11, labels updated for v17)
    -- =========================================================================
    process(clk)
        variable ip      : integer;
        variable fp      : integer;
        variable rel     : integer;
        variable raw_val : integer;
        variable diff    : integer;

        variable v_algo_total_cyc  : integer;
        variable v_init_cyc        : integer;
        variable v_update_cyc      : integer;
        variable v_forecast_cyc    : integer;
        variable v_algo_total_ns   : integer;
        variable v_init_ns         : integer;
        variable v_update_ns       : integer;
        variable v_forecast_ns     : integer;
        variable v_fitted_span_ns  : integer;
        variable v_fcast_span_ns   : integer;
        variable v_fitted_tp_whole : integer;
        variable v_fitted_tp_frac  : integer;
        variable v_fcast_tp_whole  : integer;
        variable v_fcast_tp_frac   : integer;
        variable v_e2e_tp_whole    : integer;
        variable v_e2e_tp_frac     : integer;
        variable v_total_outputs   : integer;
        variable v_init_pct        : integer;
        variable v_update_pct      : integer;
        variable v_forecast_pct    : integer;
    begin
        if rising_edge(clk) then

            if valid_in = '1' and not sv_start_snap_set then
                sv_start_snap     := cycle_count + 1;
                sv_start_snap_set := true;
                report "--- monitor: start_snap captured at cycle " &
                       integer'image(sv_start_snap)
                severity warning;
            end if;

            if fitted_valid = '1' then
                fitted_count <= fitted_count + 1;
                rel          := cycle_count + 1 - sv_start_snap;

                if not first_fitted_seen then
                    first_fitted_cycle <= rel;
                    first_fitted_seen  <= true;
                elsif not second_fitted_seen then
                    second_fitted_cycle <= rel;
                    second_fitted_seen  <= true;
                end if;

                last_fitted_cycle <= rel;
                raw_val := to_integer(fitted_out);
                decode_q1616(raw_val, ip, fp);

                if ip < 3000 or ip > 12000 then
                    report "!!! FITTED SANITY FAIL #" &
                           integer'image(fitted_count + 1) &
                           "  val=" & integer'image(ip) & "." & integer'image(fp)
                    severity warning;
                end if;

                report "FITTED  #" & integer'image(fitted_count + 1) &
                       "  rel=" & integer'image(rel) &
                       "  val=" & integer'image(ip) & "." & integer'image(fp)
                severity warning;
            end if;

            if forecast_valid = '1' then
                forecast_count <= forecast_count + 1;
                rel            := cycle_count + 1 - sv_start_snap;
                raw_val        := to_integer(forecast_out);

                if not first_forecast_seen then
                    first_forecast_cycle <= rel;
                    first_forecast_seen  <= true;
                elsif not second_forecast_seen then
                    second_forecast_cycle <= rel;
                    second_forecast_seen  <= true;
                end if;

                decode_q1616(raw_val, ip, fp);
                sv_fc_log_raw(forecast_count)   := raw_val;
                sv_fc_log_ipart(forecast_count) := ip;
                sv_fc_log_fpart(forecast_count) := fp;
                sv_fc_log_cycle(forecast_count) := rel;

                if forecast_count < 12 then
                    diff := raw_val - EXPECTED_FC(forecast_count);
                    if diff < 0 then diff := -diff; end if;
                    if diff <= FC_TOLERANCE then
                        report "FORECAST #" & integer'image(forecast_count + 1) &
                               "  rel=" & integer'image(rel) &
                               "  val=" & integer'image(ip) & "." & integer'image(fp) &
                               "  PASS (diff=" & integer'image(diff) & ")"
                        severity warning;
                    else
                        report "FORECAST #" & integer'image(forecast_count + 1) &
                               "  rel=" & integer'image(rel) &
                               "  val=" & integer'image(ip) & "." & integer'image(fp) &
                               "  FAIL diff=" & integer'image(diff) &
                               " > tol=" & integer'image(FC_TOLERANCE)
                        severity failure;
                    end if;
                end if;

                if last_forecast = '1' then
                    last_forecast_cycle <= rel;

                    v_algo_total_cyc := rel;
                    v_init_cyc       := first_fitted_cycle;
                    v_update_cyc     := last_fitted_cycle  - first_fitted_cycle;
                    v_forecast_cyc   := rel - first_forecast_cycle;
                    v_algo_total_ns  := v_algo_total_cyc * CLK_NS;
                    v_init_ns        := v_init_cyc        * CLK_NS;
                    v_update_ns      := v_update_cyc      * CLK_NS;
                    v_forecast_ns    := v_forecast_cyc    * CLK_NS;

                    if v_algo_total_cyc > 0 then
                        v_init_pct     := v_init_cyc     * 100 / v_algo_total_cyc;
                        v_update_pct   := v_update_cyc   * 100 / v_algo_total_cyc;
                        v_forecast_pct := v_forecast_cyc * 100 / v_algo_total_cyc;
                    else
                        v_init_pct := 0; v_update_pct := 0; v_forecast_pct := 0;
                    end if;

                    v_total_outputs  := fitted_count + forecast_count + 1;
                    v_fitted_span_ns := (second_fitted_cycle   - first_fitted_cycle)   * CLK_NS;
                    v_fcast_span_ns  := (second_forecast_cycle - first_forecast_cycle) * CLK_NS;

                    if v_fitted_span_ns > 0 then
                        v_fitted_tp_whole := 1000  / v_fitted_span_ns;
                        v_fitted_tp_frac  := 10000 / v_fitted_span_ns
                                             - v_fitted_tp_whole * 10;
                    else
                        v_fitted_tp_whole := 0; v_fitted_tp_frac := 0;
                    end if;

                    if v_fcast_span_ns > 0 then
                        v_fcast_tp_whole := 1000  / v_fcast_span_ns;
                        v_fcast_tp_frac  := 10000 / v_fcast_span_ns
                                            - v_fcast_tp_whole * 10;
                    else
                        v_fcast_tp_whole := 0; v_fcast_tp_frac := 0;
                    end if;

                    if v_algo_total_ns > 0 then
                        v_e2e_tp_whole := (v_total_outputs * 1000)  / v_algo_total_ns;
                        v_e2e_tp_frac  := (v_total_outputs * 10000) / v_algo_total_ns
                                          - v_e2e_tp_whole * 10;
                    else
                        v_e2e_tp_whole := 0; v_e2e_tp_frac := 0;
                    end if;

                    -- ==========================================================
                    -- Final report -- version label updated for v17
                    -- ==========================================================
                    report "===================================================" severity warning;
                    report "VERSION : v17 (beta=0,gamma=0,M=12,N=60,H=12,CLK=10ns)" severity warning;
                    report "CHANGES : DSP output pipeline registers (18 FSM states)" severity warning;
                    report "TARGET  : 100 MHz timing closure with MCP constraints" severity warning;
                    report "CONTROL : fully automatic -- valid_in only" severity warning;
                    report "REFERENCE: all cycles from first valid_in pulse" severity warning;

                    report "--- PIPELINE DEPTH vs v11 ---" severity warning;
                    report "v11: 13 states  ~355 cycles  4 cyc/sample  14ns clock" severity warning;
                    report "v17: 18 states  ~649 cycles  8 cyc/sample  10ns clock" severity warning;
                    report "extra states: INIT_LT_COMMIT MULT2 WAIT2 FINALIZE2 FORECAST_PRE" severity warning;

                    report "--- PHASE LATENCIES ---" severity warning;
                    report "Init overhead (first_vin -> first fitted_valid)" severity warning;
                    report "  " & integer'image(v_init_cyc) & " cycles  " &
                           integer'image(v_init_ns) & " ns  " &
                           "(" & integer'image(v_init_pct) & "% of total)"
                    severity warning;
                    report "Update phase (first fitted_valid -> last fitted_valid)" severity warning;
                    report "  " & integer'image(v_update_cyc) & " cycles  " &
                           integer'image(v_update_ns) & " ns  " &
                           "(" & integer'image(v_update_pct) & "% of total)"
                    severity warning;
                    report "Gap (last fitted -> first forecast)" severity warning;
                    report "  " &
                           integer'image(first_forecast_cycle - last_fitted_cycle) &
                           " cycles  " &
                           integer'image((first_forecast_cycle - last_fitted_cycle) * CLK_NS) &
                           " ns"
                    severity warning;
                    report "Forecast phase (first forecast_valid -> last forecast_valid)" severity warning;
                    report "  " & integer'image(v_forecast_cyc) & " cycles  " &
                           integer'image(v_forecast_ns) & " ns  " &
                           "(" & integer'image(v_forecast_pct) & "% of total)"
                    severity warning;

                    report "--- TOTAL ALGORITHMIC LATENCY ---" severity warning;
                    report "first_valid_in -> last forecast_valid" severity warning;
                    report "  " & integer'image(v_algo_total_cyc) & " cycles  " &
                           integer'image(v_algo_total_ns) & " ns"
                    severity warning;

                    report "--- WALL-CLOCK COMPARISON ---" severity warning;
                    report "v11 @ 71MHz: ~355 cycles x 14ns = ~4970 ns" severity warning;
                    report "v17 @ 100MHz: " &
                           integer'image(v_algo_total_cyc) & " cycles x 10ns = " &
                           integer'image(v_algo_total_ns) & " ns"
                    severity warning;

                    report "--- KEY LATENCY POINTS ---" severity warning;
                    report "first_vin -> first fitted_valid  : " &
                           integer'image(first_fitted_cycle) & " cycles  (" &
                           integer'image(first_fitted_cycle * CLK_NS) & " ns)"
                    severity warning;
                    report "first_vin -> last  fitted_valid  : " &
                           integer'image(last_fitted_cycle) & " cycles  (" &
                           integer'image(last_fitted_cycle * CLK_NS) & " ns)"
                    severity warning;
                    report "first_vin -> first forecast_valid: " &
                           integer'image(first_forecast_cycle) & " cycles  (" &
                           integer'image(first_forecast_cycle * CLK_NS) & " ns)"
                    severity warning;
                    report "first_vin -> last  forecast_valid: " &
                           integer'image(rel) & " cycles  (" &
                           integer'image(rel * CLK_NS) & " ns)"
                    severity warning;

                    report "--- ORDER CHECK ---" severity warning;
                    report "last fitted < first forecast: " &
                           integer'image(last_fitted_cycle) & " < " &
                           integer'image(first_forecast_cycle) & "  ->  " &
                           boolean'image(last_fitted_cycle < first_forecast_cycle)
                    severity warning;

                    report "--- THROUGHPUT ---" severity warning;
                    report "fitted   (update phase, inter-output period): " &
                           integer'image(v_fitted_tp_whole) & "." &
                           integer'image(v_fitted_tp_frac) & " Msamples/sec" &
                           "  (period=" & integer'image(v_fitted_span_ns) & " ns" &
                           " = " & integer'image(v_fitted_span_ns / CLK_NS) & " cycles)"
                    severity warning;
                    report "forecast (fcast phase, inter-output period): " &
                           integer'image(v_fcast_tp_whole) & "." &
                           integer'image(v_fcast_tp_frac) & " Msamples/sec" &
                           "  (period=" & integer'image(v_fcast_span_ns) & " ns" &
                           " = " & integer'image(v_fcast_span_ns / CLK_NS) & " cycles)"
                    severity warning;
                    report "end-to-end (first_vin -> last forecast): " &
                           integer'image(v_e2e_tp_whole) & "." &
                           integer'image(v_e2e_tp_frac) & " Msamples/sec" &
                           "  (" & integer'image(v_total_outputs) & " outputs" &
                           "  span=" & integer'image(v_algo_total_ns) & " ns)"
                    severity warning;

                    report "--- COUNTS ---" severity warning;
                    report "total fitted   : " & integer'image(fitted_count)       severity warning;
                    report "total forecast : " & integer'image(forecast_count + 1) severity warning;

                    report "--- FORECAST RAW LOG ---" severity warning;
                    for i in 0 to MAX_HORIZON_TB - 1 loop
                        exit when i >= forecast_count + 1;
                        report "  FC #" & integer'image(i + 1) &
                               "  rel=" & integer'image(sv_fc_log_cycle(i)) &
                               "  raw=" & integer'image(sv_fc_log_raw(i)) &
                               "  val=" & integer'image(sv_fc_log_ipart(i)) &
                               "."      & integer'image(sv_fc_log_fpart(i))
                        severity warning;
                    end loop;

                    if fitted_count = 60 and forecast_count + 1 = 12 then
                        report "TEST COMPLETE - counts OK (60 fitted, 12 forecast)"
                        severity warning;
                    else
                        report "!!! TEST FAIL  fitted=" & integer'image(fitted_count) &
                               "  forecast=" & integer'image(forecast_count + 1)
                        severity failure;
                    end if;

                    report "===================================================" severity warning;
                    done <= true;
                end if;
            end if;

            if error_out = '1' then
                report "!!! ERROR code=" &
                       integer'image(to_integer(unsigned(error_code))) &
                       "  abs=" & integer'image(abs_cycle_count)
                severity failure;
            end if;

        end if;
    end process;

    -- =========================================================================
    -- Stimulus process (protocol identical to v11 testbench)
    -- =========================================================================
    process
    begin
        rst      <= '1';
        valid_in <= '0';
        data_in  <= (others => '0');

        for i in 0 to 4 loop
            wait until rising_edge(clk);
        end loop;
        wait for CLK_PERIOD - TSETUP;
        rst <= '0';
        wait until rising_edge(clk);
        wait until rising_edge(clk);

        report "--- STEP 1: prefetch sample[0], valid_in low..." severity warning;

        -- ── STEP 1: Prefetch -- data stable, valid low ────────────────────────
        wait for CLK_PERIOD - TSETUP;
        data_in  <= PSEI_ROM(0);
        valid_in <= '0';
        wait until rising_edge(clk);

        -- ── STEP 2: First valid pulse -- S_IDLE fires, goes to S_COLLECT ─────
        wait for CLK_PERIOD - TSETUP;
        data_in  <= PSEI_ROM(0);
        valid_in <= '1';
        report "--- STEP 2: first valid_in asserted (S_IDLE -> S_COLLECT)" severity warning;
        wait until rising_edge(clk);

        -- ── STEP 3: Hold ROM(0) -- S_COLLECT writes it to data_buf(0) ────────
        wait for CLK_PERIOD - TSETUP;
        data_in  <= PSEI_ROM(0);
        valid_in <= '1';
        report "--- STEP 3: holding ROM(0) -- S_COLLECT writes sample[0]" severity warning;
        wait until rising_edge(clk);

        -- ── STEP 4: Stream samples 1..59 ──────────────────────────────────────
        report "--- STEP 4: streaming samples 1..59..." severity warning;
        for i in 1 to ROM_SIZE - 1 loop
            wait for CLK_PERIOD - TSETUP;
            data_in  <= PSEI_ROM(i);
            valid_in <= '1';
            wait until rising_edge(clk);
        end loop;

        -- ── Drop valid_in ──────────────────────────────────────────────────────
        wait for CLK_PERIOD - TSETUP;
        valid_in <= '0';
        data_in  <= (others => '0');
        report "--- stream complete (60 samples). Core self-drives." severity warning;

        -- ── Wait for completion ────────────────────────────────────────────────
        wait until done = true;
        wait until rising_edge(clk);
        wait until rising_edge(clk);

        std.env.stop;
        wait;
    end process;

end architecture sim;