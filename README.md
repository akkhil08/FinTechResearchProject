# GA-Optimised Holt-Winters Forecaster — FPGA Hardware Accelerator

A complete end-to-end hardware-software co-design project implementing
Genetic Algorithm-optimised Holt-Winters triple exponential smoothing,
accelerated on the PL (Programmable Logic) fabric of a Xilinx Zynq-7000
ZedBoard. Parameters are optimised offline in C++ and hardcoded into a
VHDL forecasting core — no ARM PS involvement at runtime.

---

## Results

| Metric                     | Value                                  |
|----------------------------|----------------------------------------|
| Target device              | XC7Z020CLG484-1 (Zynq-7000 ZedBoard)  |
| Implementation target      | PL fabric only (no PS at runtime)      |
| Clock frequency            | 100 MHz (10 ns period)                 |
| Worst Negative Slack (WNS) | +0.537 ns                              |
| End-to-end latency         | 651 cycles / 6.51 µs                   |
| CPU average latency (1000 runs) | 18,434 ns / 18.43 µs              |
| FPGA speedup vs CPU        | 2.83× faster — 11,934 ns saved per forecast |
| FSM states                 | 18                                     |
| Parameter format           | Q2.30 fixed-point (α, β, γ)           |
| Data format                | Q16.16 fixed-point (level, trend, seasonal, forecast) |
| Optimisation versions      | 13 iterative VHDL revisions            |

### Resource utilisation — XC7Z020CLG484-1

| Resource               | Used  | Available | Utilisation |
|------------------------|-------|-----------|-------------|
| Slice LUTs (Logic)     | 2,141 | 53,200    | 4.0%        |
| LUT as Distributed RAM | 220   | 17,400    | 1.3%        |
| Slice Registers (FFs)  | 1,702 | 106,400   | 1.6%        |
| DSP48E1                | 56    | 220       | 25.5%       |
| Block RAM Tile         | 0.5   | 140       | 0.4%        |

### GA optimisation results

| Model          | Training MSE | Selected |
|----------------|-------------|----------|
| Additive       | 108,739     | ✓        |
| Multiplicative | 113,583     |          |

GA-converged parameters: **α = 0.9516, β = 0, γ = 0**

The zero values of β and γ indicate the GA found negligible trend and seasonal
components over the 2017–2021 training window, effectively reducing the model
to a first-order exponential smoother with fixed initial seasonal structure.

### Verification — C++ vs VHDL (PSEI 2022)

| Month | Actual  | C++ Forecast | VHDL Forecast | \|C++−VHDL\| |
|-------|---------|-------------|---------------|--------------|
| Jan   | 7361.65 | 7062.108    | 7062.111      | 0.003        |
| Feb   | 7311.01 | 6881.042    | 6881.051      | 0.009        |
| Mar   | 7203.47 | 6467.992    | 6468.010      | 0.018        |
| Apr   | 6731.25 | 6560.262    | 6560.282      | 0.020        |
| May   | 6774.68 | 6600.396    | 6600.418      | 0.022        |
| Jun   | 6155.43 | 6662.022    | 6662.046      | 0.024        |
| Jul   | 6315.93 | 6606.410    | 6606.438      | 0.028        |
| Aug   | 6583.65 | 6712.934    | 6712.964      | 0.030        |
| Sep   | 5741.07 | 6601.826    | 6601.860      | 0.034        |
| Oct   | 6153.43 | 6751.932    | 6751.967      | 0.035        |
| Nov   | 6780.78 | 6837.012    | 6837.049      | 0.037        |
| Dec   | 6566.39 | 6973.503    | 6973.540      | 0.037        |
| **MAE** |       | **388.27**  | **388.28**    | **max 0.037** |

VHDL outputs are numerically consistent with the C++ baseline within
Q16.16 fixed-point quantisation bounds across all forecast horizons.

![C++ vs VHDL Forecast Comparison](results/vhdl_vs_cpp_forecast.png)

---

## Project Overview

### Problem
Real-time financial embedded systems require deterministic, low-latency
forecasting. CPU-based Holt-Winters operates in the millisecond range —
unsuitable for time-critical applications. FPGA deployment on the PL fabric
reduces this to the microsecond range with fully deterministic behaviour.

### Solution

**Stage 1 — C++ software pipeline**

A self-contained C++ programme parses PSEI historical closing price data,
applies a Genetic Algorithm to find optimal Holt-Winters smoothing parameters
(α, β, γ ∈ [0,1]), evaluates both additive and multiplicative models, selects
the best on training MSE, and exports verified forecast outputs to CSV as the
numerical baseline.

GA configuration: population 40, generations 80, tournament selection (k=3),
arithmetic crossover (rate 0.8), Gaussian mutation (σ=0.05, rate 0.2).

**Stage 2 — VHDL hardware accelerator (PL fabric)**

The GA-optimised parameters are hardcoded into a VHDL design targeting the
XC7Z020CLG484-1 PL fabric. Two modules:

- `hw_top` — stores 60 PSEI training values as compile-time constants in
  Block RAM (RAM64M) and streams them into the forecasting core via valid_in/data_in
  handshake which makes copies according to number of address variables in the core. 
- `holt_winters_q2_30` — implements the full Holt-Winters additive algorithm
  as an 18-state FSM using Q2.30 for smoothing parameters and Q16.16 for data
  values; all multiplications use DSP48E1 primitives

**Stage 3 — Numerical verification**

VHDL simulation outputs are compared against the C++ baseline on 2022 PSEI
data. Maximum deviation: 0.037 index points, confirming implementation
correctness within Q16.16 fixed-point quantisation bounds.

---

## Architecture

```
┌──────────────────────────────────────┐
│       GA Optimiser (C++ offline)     │
│  Population: 40 · Generations: 80    │
│  Selection: Tournament (k=3)         │
│  Crossover: Arithmetic (rate=0.8)    │
│  Mutation:  Gaussian (σ=0.05)        │
│  Output: α=0.9516, β=0, γ=0         │
└──────────────┬───────────────────────┘
               │ hardcoded parameters
               ▼
┌──────────────────────────────────────┐   XC7Z020CLG484-1
│         PL Fabric (VHDL)            │
│                                      │
│  hw_top                              │
│  └── holt_winters_q2_30             │
│      ├── 18-state FSM datapath       │
│      ├── Q2.30 params (α, β, γ)     │
│      ├── Q16.16 data arithmetic      │
│      ├── DSP48E1 multipliers (×56)   │
│      └── ILA debug core (v2 only)    │
│                                      │
│  PS (ARM Cortex-A9): NOT used        │
└──────────────┬───────────────────────┘
               │ 651 cycles / 6.51 µs
               ▼
┌──────────────────────────────────────┐
│     Verification (Python/Jupyter)    │
│     Max deviation: 0.037 pts         │
│     MAE: 388.27 (C++) / 388.28 (HW) │
└──────────────────────────────────────┘
```

---

## Hardware Optimisation Journey

13 VHDL versions were required to progress from 48.78 MHz to 100 MHz closure.

| Ver. | Key Change | LUTs (Logic/Total) | FFs (Logic/Total) | DSPs | WNS (ns) | Clock | Cycles | Latency (ns) |
|------|-----------|-------------------|------------------|------|----------|-------|--------|--------------|
| v1   | Hardcoded baseline | 3260/3180 | 2323/2323 | 26 | +0.249 | 48.78 MHz | — | 50,505 |
| ILA  | On-board ILA verification | 3671/4024 | 2100/2981 | 34 | −31.882 | — | — | — |
| v2   | Generics + 8 RAM replicas | 5442/5412 | 2754/2853 | 88 | −28.106 | 71.42 MHz | 217 | 3,038 |
| v3   | Odd/even state split | 5404/5365 | 2745/2844 | 88 | −26.973 | 71.42 MHz | 223 | 3,122 |
| v4   | Ping-pong L/T registers | 5626/5594 | 2763/2862 | 88 | −27.513 | 71.42 MHz | 223 | 3,122 |
| v5   | Real-time calc, 5 replicas | 4492/4467 | 2713/2812 | 86 | −26.231 | 71.42 MHz | 283 | 3,962 |
| v6   | Pre-computed register caching | 5623/5572 | 2734/2788 | 138 | −17.316 | 71.42 MHz | 347 | 4,858 |
| v8   | Pipeline + DSP targeting | 4025/3939 | 3031/3034 | 80 | +0.335 | 71.42 MHz | 546 | 7,644 |
| v9   | Pipeline latency recovery | 4248/4160 | 3026/3029 | 80 | +0.293 | 71.42 MHz | 487 | 6,818 |
| v10  | Single address, replicas removed | 2966/2860 | 2193/2196 | 80 | +0.198 | 71.42 MHz | 403 | 5,642 |
| v11  | Integer constants, DSP reduction | 2119/2029 | 1396/1399 | 56 | +0.449 | 71.42 MHz | 403 | 5,642 |
| v12  | State and variable cleanup | 2081/2017 | 1374/1389 | 56 | +0.287 | 71.42 MHz | 458 | 6,412 |
| **v13** | **100 MHz closure** | **2429/2361** | **1702/1702** | **56** | **+0.537** | **100 MHz** | **651** | **6,510** |

The `hardware/archive/` folder preserves all 13 versions.

Key engineering challenges:

- **v2–v5 RAM replica problem** — generalising the design forced Vivado to replicate
  `data_buf` for each unique read address signal. Confirmed via TCL console: 55 RAM64M
  instances across 5 replicas. Resolved in v10 by consolidating all reads to a single
  registered address.

- **v8 first closure at 71 MHz** — splitting the update recurrence into one-multiply-per-state
  pipeline stages with explicit `use_dsp` attributes achieved WNS +0.335 ns.

- **v13 100 MHz closure** — DSP48E1 outputs on XC7Z020-1 consume ~8.5 ns internally,
  leaving <1.5 ns for routing. Resolved by adding a raw-register hop state (S_UPDATE_WAIT2)
  adjacent to the DSP primitives, plus multicycle path constraints (2-cycle setup / 1-cycle
  hold) in the XDC file.

---

## Repository Structure

```
├── README.md
│
├── software/
│   ├── HWlogideawithGA.cpp                # GA-optimised Holt-Winters C++ pipeline
│   └── holtwinter_timing_measurement.cpp  # Software timing benchmark
│
├── hardware/
│   ├── hw_top.vhd                         # Top-level: ROM streaming + handshake
│   ├── hw_version13.vhd                   # Final forecasting core (holt_winters_q2_30)
│   ├── hw_version13_tb.vhd                # Testbench — timing-verified simulation
│   └── archive/                           # All 13 iterative versions (v1–v12)
│
├── verification/
│   ├── hw_ila_idea_check.vhd              # ILA signal capture (used in v2)
│   └── hw_ila_idea_check_top.vhd         # ILA top-level integration
│
└── results/
    ├── VHDLC++.ipynb                      # Jupyter notebook — verification analysis
    ├── c++output                          # Raw C++ software output
    ├── combined_forecast_data.csv         # Merged C++ and VHDL forecast data
    ├── cpu_vs_fpga_timing.txt            # CPU vs FPGA latency comparison (1000 runs)
    ├── vhdl_simulation_output.txt         # XSim log — 651 cycles, TEST COMPLETE
    └── vhdl_vs_cpp_forecast.png           # C++ vs VHDL comparison plot
```

---

## Dataset

**Philippine Stock Exchange Index (PSEi)**
- Source: [Investing.com PSEI Historical Data](https://www.investing.com/indices/psei-composite-historical-data)
- Training: 2017–2021 monthly closing values (N=60, M=12, 5 seasons)
- Validation: 2022 monthly closing values (12-month forecast horizon)
- Live records: [Google Sheets](https://docs.google.com/spreadsheets/d/1P3ti_ZSYGuJDgwotsEgAts41fRzfiONRrHpxk0O3xAo/edit?usp=sharing)

Training window restricted to 2017–2021 (vs full 2010–2021) to avoid
structural breaks from the 2015 correction and COVID-19 volatility of 2020,
reducing MAE from ~1600 to ~388 index points.

---

## How to Build and Run

### Requirements
- Xilinx Vivado 2023.x (synthesis, implementation, XSim)
- Xilinx ZedBoard (XC7Z020CLG484-1)
- GCC with C++11 (software baseline)
- Python 3.x — pandas, numpy, matplotlib (verification plots)

### Software baseline (C++)

```bash
cd software
g++ -O2 -std=c++11 -o hw_ga HWlogideawithGA.cpp -lm
./hw_ga
# Outputs: forecast_2017_2021_vs_2022_add_vs_mul_log_vs_orig_GA.csv
```

### VHDL simulation (Vivado XSim)

1. Open Vivado → Create project
2. Add `hardware/hw_version13.vhd` (source) and `hardware/hw_version13_tb.vhd` (sim)
3. Run Behavioural Simulation
4. Verify against `results/vhdl_simulation_output.txt` — expect TEST COMPLETE, 651 cycles

### On-hardware deployment

1. Add all files in `hardware/` to Vivado project
2. Apply XDC constraints (multicycle path for DSP outputs, pin assignments)
3. Synthesis → Implementation → Generate Bitstream
4. Program ZedBoard PL via JTAG
5. Use ILA cores in `verification/` for on-hardware signal capture

### Verification plot

```bash
cd results
jupyter notebook VHDLC++.ipynb
# Reproduces vhdl_vs_cpp_forecast.png
```

---

## Fixed-Point Arithmetic

Two Q-format representations are used:

| Format | Used for | Precision | Range |
|--------|----------|-----------|-------|
| Q2.30  | α, β, γ smoothing parameters | 2⁻³⁰ ≈ 9.3×10⁻¹⁰ | [0, 1] |
| Q16.16 | Data values, level, trend, seasonal, forecast | 2⁻¹⁶ ≈ 1.5×10⁻⁵ | ±32,767 |

Multiplying Q2.30 × Q16.16 yields a Q18.46 intermediate result, right-shifted
by 30 bits to recover Q16.16. Parameters were initially Q16.16 but switched to
Q2.30 in later versions, which improved numerical consistency between C++ and
VHDL outputs.

---

## Project Context

- **Degree**: M.Sc. Electrotechnik (Electronics Design Technology)
- **University**: Universität Siegen, Germany
- **Supervisor**: M.Sc. Aravinda Lasya Indukuri
- **Chair**: Embedded Systems, University of Siegen
- **Submitted**: June 3rd, 2026
- **Scope**: Research project as part of 9 credits in course curriculum
---

## Author

**Akkhilesh Raghuram**

 [linkedin.com/in/akkhilesh-raghuram](https://linkedin.com/in/akkhilesh-raghuram)
