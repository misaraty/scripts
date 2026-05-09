## [中文版本](https://www.misaraty.com/2026-05-03_juxaid/)

## JUXAID

`PYXAID`, developed by `Oleg V. Prezhdo` and `Alexey V. Akimov`, is a well-established nonadiabatic molecular dynamics software widely used for excited-state simulations in condensed matter systems. Subsequently, `Libra` developed by `Alexey V. Akimov` further extended the methodological framework. In addition, `Hefei-NAMD`, `NEXMD`, `SHARC`, and `Newton-X` have also been widely applied in excited-state dynamics studies and demonstrated strong performance.

Based on `PYXAID`, this project develops `JUXAID`, a lightweight reimplementation in `Julia`. The code is reconstructed following the original program logic and adopts a concise single-file structure. While retaining the core ab initio electron–nuclear coupling formalism, it leverages Julia’s high-performance numerical capabilities and cross-platform execution to achieve efficient nonadiabatic molecular dynamics simulations. `JUXAID` further incorporates improved surface hopping algorithms (e.g., the `SDM` method) and a more flexible modular design, enhancing its applicability and extensibility for complex systems. It supports both `Linux` and `Windows` platforms and is suitable for large-scale simulations as well as future integration with machine learning methods.

## Usage

`julia namd.jl`

## Benchmark

Performance benchmarking was conducted for `PYXAID`, `MAXAID`, and `JUXAID` (different versions) on the following platforms:

- Server platform: `Intel® Xeon® Gold 6444Y`, 256 GB RAM, `CentOS`

- Desktop platform: `Intel(R) Core(TM) Ultra 9 285H`, 32 GB RAM, `Windows 11` 

### Relative Speed Comparison

| Software Version | Relative Speed (`PYXAID` = 1) | Performance Description |
|---|---|---|
| **`PYXAID`** | **1×** | **Baseline reference** |
| `MAXAID` | 1/15× | Performance decreases significantly as the number of electronic states increases |
| `JUXAID (v10)` | 1× | Comparable to `PYXAID` |
| `JUXAID (v12)` | 6× | Optimized |
| `JUXAID (v14)` | 13× | Further optimized |
| **`JUXAID (v15)`** | **19×** | **Partial non-essential outputs disabled** |
| `PXAID (v10)` | 1/36× | Initial Python reimplementation with relatively low performance |
| `PXAID (v12)` | 1/15× | Optimized |
| `PXAID (v13)` | 11× | Optimized |
| `PXAID (v14)` | 34× | Further optimized |
| **`PXAID (v15)`** | **154×** | **Partial non-essential outputs disabled** |

## Changelog

### Major Improvements of v15 Compared to v14

#### 1. Added Output File Save Switches

- v14 saved all output files by default

- v15 introduces four output control parameters:

  - `SAVE_DECOHERENCE_RATES`
  
  - `SAVE_ICOND_FILES`
  
  - `SAVE_ME_ENERGIES`
  
  - `SAVE_ME_POP`

- Users can flexibly choose whether to save specific output files

#### 2. Reduced File IO Overhead

- v14 continuously wrote multiple text output files

- v15 allows disabling non-essential outputs

- Significantly reduces file writing time in large-scale `icond` and long-time trajectory simulations

### Major Improvements of v14 Compared to v12

#### 1. Introduced Multithreaded Parallel Computing (`Threads`)

- v12 mainly used single-thread execution

- v14 introduced `Base.Threads`

- Parallel acceleration was implemented for decoherence calculations and multi-trajectory dynamics simulations

- Significantly improved efficiency for large-scale trajectory calculations

#### 2. Fully “De-objectified” Core Electronic Dynamics

- v12 propagated dynamics using `ElectronicStructure` objects

- v14 rewrote the core propagation using direct low-level array operations on:

  - `C`
  
  - `A`
  
  - `g`
  
  - `tau_m`
  
  - `t_m`

- Reduced object access and dynamic dispatch overhead

#### 3. Extensive Function `inline` Optimization

- v12 contained many ordinary function calls

- v14 applied `@inline` to core small functions

- Reduced function call overhead and improved hotspot loop performance

#### 4. Batch Storage of `Hamiltonian` Using 3D Arrays

- v12 stored `Hme_batch` using `Vector{Matrix}`

- v14 replaced it with:

  - `Array{ComplexF64,3}`

- Improved memory locality and cache efficiency

#### 5. Further Compression of Orbital Mapping

- v12 used nested `Vector` structures

- v14 introduced:

  - `pack_diag_orbital_map`

- Compressed orbital mappings into fixed-size matrices

- Reduced dynamic memory access overhead

#### 6. Lower-Level Propagation Core Implementation

- v12 used higher-level object interfaces

- v14 split the implementation into multiple `_core!` functions:

  - `propagate_coefficients_core!`
  
  - `update_populations_core!`
  
  - `hop_core!`
  
  - `sdm_decoherence_core!`

- More suitable for compiler optimization and parallel execution

#### 7. Further FFT Optimization

- v12 used:

  - `fft`
  
  - `ifft`

- v14 replaced them with:

  - `rfft`
  
  - `irfft`

  in real-valued decoherence autocorrelation calculations

- Reduced FFT computational cost and memory usage by exploiting real-valued energy fluctuation sequences

#### 8. Optimized Large-Scale Output Writing

- v12 performed frequent direct file writes

- v14 introduced:

  - `IOBuffer`

- Reduced performance loss caused by大量 small-scale IO operations

#### 9. Thread-Localized Multi-Trajectory Simulation

- v12 used globally shared population arrays for all trajectories

- v14 introduced per-thread local arrays:

  - `sh_tls`
  
  - `se_tls`

- Reduced thread contention and synchronization overhead

#### 10. Further Optimization of Data Loading Toward Contiguous Memory

- v12 stored `H_batch` using `Vector{Matrix}`

- v14 replaced it with:

  - `H_array`

- Better suited for batched linear algebra operations and cache optimization

### Major Improvements from v10 to v12

#### 1. FFT-based Autocorrelation

- v10 computes autocorrelation using nested loops (O(N²))

- v12 uses FFT-based implementation (O(N log N))

- Significantly accelerates decoherence calculations for long trajectories

#### 2. Parameter Struct (Remove Globals)

- v10 relies heavily on global variables (`hbar`, `kb`, `dt`, etc.)

- v12 introduces `NAMDParams` struct

- Improves type stability and compiler optimization

#### 3. Precomputed Many-Electron Mapping

- v10 repeatedly computes `delta_states` and `ext2int`

- v12 precomputes:

  - `transition_map`
  
  - `diag_orbital_map`
  
- Reduces redundant computations

#### 4. Avoid Full Time-Series Object Storage

- v10 stores full `oe_es` and `me_es` arrays (memory intensive)

- v12 uses:

  - `Hme_batch` (lightweight)
  
  - Single `ElectronicStructure` object
  
- Greatly reduces memory usage

#### 5. Modular Hamiltonian Construction

- v10 builds Hamiltonian inside `main`

- v12 refactors into:

  - `build_spin_orbital_H!`
  
  - `build_Hme_batch`
  
  - `write_me_energies`
  
- Improves readability and maintainability

#### 6. Reduced Memory Allocation

- v10 uses broadcasting and matrix multiplication (temporary arrays)

- v12 uses explicit loops

- Reduces `GC` overhead and improves performance

#### 7. Optimized Hop Algorithm

- v10 uses `vec` + `cumsum`

- v12 uses `view` + incremental accumulation

- Avoids unnecessary allocations

#### 8. Parallelization-Friendly Design

- v12 uses flattened data structures (`Hme_batch` + single object)

- Easier to extend for:

  - multithreading
  
  - distributed computing

## Citation

Original `PYXAID` references:

* Akimov A V, Prezhdo O V. The PYXAID program for non-adiabatic molecular dynamics in condensed matter systems. Journal of Chemical Theory and Computation, 2013, 9(11): 4959–4972.

* Akimov A V, Prezhdo O V. Advanced capabilities of the PYXAID program: integration schemes, decoherence effects, multiexcitonic states, and field–matter interaction. Journal of Chemical Theory and Computation, 2014, 10(2): 789–804.

This work:

To be added after the paper is officially published.