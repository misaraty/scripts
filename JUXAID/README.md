## [中文版本](https://www.misaraty.com/2026-04-30_juxaid/)

## JUXAID

`PYXAID`, developed by `Oleg V. Prezhdo` and `Alexey V. Akimov`, is a well-established nonadiabatic molecular dynamics software widely used for excited-state simulations in condensed matter systems. Subsequently, `Libra` developed by `Alexey V. Akimov` further extended the methodological framework. In addition, `Hefei-NAMD`, `NEXMD`, `SHARC`, and `Newton-X` have also been widely applied in excited-state dynamics studies and demonstrated strong performance.

Based on `PYXAID`, this project develops `JUXAID`, a lightweight reimplementation in `Julia`. The code is reconstructed following the original program logic and adopts a concise single-file structure. While retaining the core ab initio electron–nuclear coupling formalism, it leverages Julia’s high-performance numerical capabilities and cross-platform execution to achieve efficient nonadiabatic molecular dynamics simulations. `JUXAID` further incorporates improved surface hopping algorithms (e.g., the `SDM` method) and a more flexible modular design, enhancing its applicability and extensibility for complex systems. It supports both `Linux` and `Windows` platforms and is suitable for large-scale simulations as well as future integration with machine learning methods.

## Usage

`julia namd.jl`

## Benchmark

Performance benchmarking was conducted for `PYXAID`, `MAXAID`, and `JUXAID` (different versions) on the following platforms:

- Server platform: `Intel® Xeon® Gold 6444Y`, 256 GB RAM, `CentOS`

- Desktop platform: `Intel(R) Core(TM) Ultra 9 285H`, 32 GB RAM, `Windows 11` 

### Relative Runtime Comparison

| Software        | Relative Runtime (`PYXAID`=1) | Remarks |
|----------------|-----------------------------|--------|
| `PYXAID`         | 1.0                         | Baseline |
| `JUXAID (v9)`    | ~1.0                        | Comparable to `PYXAID` |
| `MAXAID`         | 10 ~ 20                     | Significant slowdown with increasing states |
| `JUXAID (v11)`   | ~0.17                       | ~6× speedup over `PYXAID` |

## Changelog

### Major Improvements from v9 to v11

#### 1. FFT-based Autocorrelation

- v9 computes autocorrelation using nested loops (O(N²))

- v11 uses FFT-based implementation (O(N log N))

- Significantly accelerates decoherence calculations for long trajectories

#### 2. Parameter Struct (Remove Globals)

- v9 relies heavily on global variables (`hbar`, `kb`, `dt`, etc.)

- v11 introduces `NAMDParams` struct

- Improves type stability and compiler optimization

#### 3. Precomputed Many-Electron Mapping

- v9 repeatedly computes `delta_states` and `ext2int`

- v11 precomputes:

  - `transition_map`
  
  - `diag_orbital_map`
  
- Reduces redundant computations

#### 4. Avoid Full Time-Series Object Storage

- v9 stores full `oe_es` and `me_es` arrays (memory intensive)

- v11 uses:

  - `Hme_batch` (lightweight)
  
  - Single `ElectronicStructure` object
  
- Greatly reduces memory usage

#### 5. Modular Hamiltonian Construction

- v9 builds Hamiltonian inside `main`

- v11 refactors into:

  - `build_spin_orbital_H!`
  
  - `build_Hme_batch`
  
  - `write_me_energies`
  
- Improves readability and maintainability

#### 6. Reduced Memory Allocation

- v9 uses broadcasting and matrix multiplication (temporary arrays)

- v11 uses explicit loops

- Reduces `GC` overhead and improves performance

#### 7. Optimized Hop Algorithm

- v9 uses `vec` + `cumsum`

- v11 uses `view` + incremental accumulation

- Avoids unnecessary allocations

#### 8. Parallelization-Friendly Design

- v11 uses flattened data structures (`Hme_batch` + single object)

- Easier to extend for:

  - multithreading
  
  - distributed computing

## Citation

Original `PYXAID` references:

* Akimov A V, Prezhdo O V. The PYXAID program for non-adiabatic molecular dynamics in condensed matter systems. Journal of Chemical Theory and Computation, 2013, 9(11): 4959–4972.

* Akimov A V, Prezhdo O V. Advanced capabilities of the PYXAID program: integration schemes, decoherence effects, multiexcitonic states, and field–matter interaction. Journal of Chemical Theory and Computation, 2014, 10(2): 789–804.

This work:

To be added after the paper is officially published.