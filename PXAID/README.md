## [中文版本](https://www.misaraty.com/2026-05-04_pxaid/)

## PXAID

`PYXAID`, developed by `Oleg V. Prezhdo` and `Alexey V. Akimov`, is a well-established nonadiabatic molecular dynamics software widely used for excited-state simulations in condensed matter systems. Subsequently, `Libra` developed by `Alexey V. Akimov` further extended the methodological framework. In addition, `Hefei-NAMD`, `NEXMD`, `SHARC`, and `Newton-X` have also been widely applied in excited-state dynamics studies and demonstrated strong performance.

Based on `PYXAID`, this project develops `PXAID`, a lightweight reimplementation written in the `Python` language. The implementation reconstructs the original program logic while adopting a more modern `Python` numerical computing framework. By preserving the core theoretical framework of *ab initio* electron–nuclear coupled dynamics, `PXAID` fully leverages high-performance scientific computing libraries such as `NumPy`, `SciPy`, and `Numba` to achieve efficient nonadiabatic molecular dynamics simulations. Key computational modules, including electronic propagation, surface hopping, and decoherence treatments, are accelerated using `Numba JIT` compilation, while `FFT`-based approaches are employed to optimize autocorrelation and decoherence-rate calculations, significantly improving computational efficiency for long-time trajectories. In addition, `PXAID` incorporates improved surface hopping algorithms (e.g., the `SDM` method) together with a more flexible modular design, enhancing its applicability and extensibility for complex systems. Furthermore, `PXAID` supports cross-platform execution on environments such as `Linux` and `Windows`, making it suitable for large-scale simulations and future integration with machine learning approaches.

## Usage

`python namd.py`

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

## Update Log

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

### Major Improvements of v14 Compared to v13

#### 1. Introduced Multithreaded Parallel Computing (`Numba parallel`)

- v13 mainly used single-thread trajectory propagation

- v14 introduced:

  - `parallel=True`
  
  - `prange`

  in the core propagation functions

- Enabled parallel computation of multiple surface hopping trajectories

- Significantly improved efficiency for large-scale trajectory simulations

#### 2. Introduced Thread-Local Storage Arrays (TLS)

- v13 shared global arrays for all trajectories:

  - `sh_pops`
  
  - `se_pops`

- v14 introduced per-thread local arrays:

  - `sh_tls`
  
  - `se_tls`

- Reduced multithreaded write conflicts and synchronization overhead

#### 3. Static Multi-Trajectory Task Partitioning

- v13 sequentially iterated over:

  - `for itraj in range(num_sh_traj)`

- v14 statically partitioned trajectory ranges across threads using:

  - `lo`
  
  - `hi`

- Improved thread load balancing and cache locality

#### 4. Independent Random Number Streams for Multithreading

- v13 shared a global random sequence across trajectories

- v14 used:

  - `seed + 1000003 * tid`

- Established independent random streams for different threads

- Reduced correlations in parallel random number generation

#### 5. Fully Parallelized Core NAMD Loop

- v13 only JIT-compiled core numerical kernels

- v14 integrated:

  - surface hopping
  
  - electronic propagation
  
  - decoherence
  
  - population accumulation

  into a fully parallel JIT core

- Further reduced Python-level scheduling overhead

#### 6. Reduced Shared Array Contention

- v13 directly updated global population arrays

- v14 used thread-local accumulation followed by reduction

- Improved multicore scalability

#### 7. Further Integration of Parallelization and JIT

- v13 used:

  - `@njit(cache=True, fastmath=True)`

- v14 changed the core propagation function to:

  - `@njit(fastmath=True, parallel=True)`

- Further improved hotspot loop execution efficiency

#### 8. Preserved Flattened Array Structures in Multithreaded Execution

- v13 already used:

  - `np.ndarray[:,:, :]`

  for batch storage of `Hamiltonian`

- v14 preserved contiguous memory layouts under parallel execution

- Improved cache efficiency and thread data access performance

#### 9. Further Reduced Python Interpreter Involvement

- v13 still retained partial Python-level trajectory control

- v14 fully moved multi-trajectory propagation into the parallel JIT core

- Reduced interpreter and function scheduling overhead

#### 10. Significantly Improved Overall Performance

- v13 achieved approximately:

  - `11×`

  relative speed compared to `PYXAID`

- v14 further improved performance under multithreaded execution

- Achieved higher performance than `JUXAID v14` in the current benchmark systems

### Major Improvements of v13 Compared to v12

#### 1. Introduced `Numba JIT` Acceleration for Core Numerical Computations

- v12 mainly relied on Python interpreted execution

- v13 introduced:

  - `@njit(cache=True, fastmath=True)`

  for core hotspot functions

- Significantly reduced Python interpreter overhead

#### 2. Fully JIT-Compiled Core Electronic Dynamics

- v12 propagated dynamics using ordinary Python functions

- v13 compiled core modules including:

  - electronic propagation
  
  - surface hopping
  
  - decoherence
  
  - population updates

  directly into machine code

- Greatly improved long-time dynamics simulation efficiency

#### 3. Further FFT Optimization

- v12 used:

  - `np.fft.fft`
  
  - `np.fft.ifft`

- v13 replaced them with:

  - `scipy.fft.rfft`
  
  - `scipy.fft.irfft`

  in real-valued decoherence autocorrelation calculations

- Reduced FFT computational cost and memory usage using real-valued energy fluctuation sequences

#### 4. Introduced `next_fast_len` Optimization for FFT Length

- v12 used fixed power-of-two zero padding

- v13 introduced:

  - `next_fast_len`

- Automatically selected more efficient FFT sizes

#### 5. Batch Storage of `Hamiltonian` Using 3D Arrays

- v12 used `list[np.ndarray]`

- v13 replaced it with:

  - `np.ndarray[:,:, :]`

- Improved memory locality and cache efficiency

#### 6. “De-objectified” Core Propagation Functions

- v12 propagated dynamics using `ElectronicStructure` objects

- v13 directly manipulated low-level arrays such as:

  - `C`
  
  - `A`
  
  - `g`
  
  - `tau_m`
  
  - `t_m`

  within JIT kernels

- Reduced object access and dynamic dispatch overhead

#### 7. Reduced Temporary Arrays and Python-Level Loops

- v12 still contained many Python-level loops

- v13 extensively adopted explicit in-place array updates inside JIT kernels

- Significantly reduced memory allocation and GC overhead

#### 8. Further JIT Optimization of `Hamiltonian` Construction

- v12 constructed `Hme_batch` using ordinary Python functions

- v13 introduced:

  - `_build_spin_orbital_H_jit`
  
  - `_build_Hme_batch_jit`

- Accelerated multielectronic Hamiltonian construction

#### 9. JIT-Compiled Random Surface Hopping Processes

- v12 used Python `random`

- v13 used:

  - `np.random`

  integrated directly into JIT kernels

- Improved surface hopping simulation efficiency

#### 10. Significantly Improved Overall Performance

- v12 achieved approximately:

  - `0.067×`

  relative speed compared to `PYXAID`

- v13 improved performance to:

  - `11×`

  relative speed compared to `PYXAID`

- Achieved order-of-magnitude performance improvement

### Major Improvements of v12 Compared to v10

#### 1. FFT Acceleration for Autocorrelation Functions

- v10 used double loops to compute autocorrelation functions (O(N²))

- v12 implemented mathematically equivalent FFT-based autocorrelation calculations (O(N log N))

- Significantly improved `decoherence` calculation efficiency for long trajectories

#### 2. Structured Parameter Management (Removal of Large Numbers of Global Variables)

- v10 used大量 global variables (`hbar`, `kb`, `dt`, etc.)

- v12 introduced the `NAMDParams` data class for unified parameter management

- Improved code maintainability and function interface consistency

#### 3. Precomputation of Multi-Electronic-State Mappings

- v10 repeatedly computed electronic state mappings at every timestep

- v12 precomputed:

  - `transition_map`
  
  - `diag_orbital_map`

- Significantly reduced repeated computational overhead

#### 4. Modularized `Hamiltonian` Construction

- v10 constructed `Hamiltonian` directly inside the main program

- v12 split the implementation into multiple independent functions

- Improved readability and extensibility

#### 5. Avoided Frequent Object Copies

- v10 propagated multiple `ElectronicStructure` object lists

- v12 switched to single-object incremental propagation

- Reduced memory usage and object management overhead

#### 6. Optimized Array Operations

- v10 contained many Python-level loops and temporary array allocations

- v12 adopted more `NumPy` vectorized operations

- Improved matrix operation efficiency

#### 7. Optimized `hop` Algorithm

- v10 used:

  - `copy`
  
  - `cumsum`
  
  - `where`

- v12 switched to online cumulative probability accumulation

- Reduced temporary array allocations

#### 8. Better Foundation for Subsequent High-Performance Optimization

- v12 further standardized data structures

- Provided the basis for later:

  - `Numba JIT`
  
  - parallel computing
  
  - FFT optimization

## Citation

Original `PYXAID` references:

* Akimov A V, Prezhdo O V. The PYXAID program for non-adiabatic molecular dynamics in condensed matter systems. Journal of Chemical Theory and Computation, 2013, 9(11): 4959–4972.

* Akimov A V, Prezhdo O V. Advanced capabilities of the PYXAID program: integration schemes, decoherence effects, multiexcitonic states, and field–matter interaction. Journal of Chemical Theory and Computation, 2014, 10(2): 789–804.

This work:

To be added after the paper is officially published.