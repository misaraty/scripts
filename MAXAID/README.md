## [中文版本](https://www.misaraty.com/2026-04-30_maxaid/)

## MAXAID

`PYXAID`, developed by `Oleg V. Prezhdo` and `Alexey V. Akimov`, is a well-established nonadiabatic molecular dynamics software widely used in excited-state simulations of condensed matter systems. `Libra` further extends this framework. In addition, `Hefei-NAMD`, `NEXMD`, `SHARC`, and `Newton-X` are also widely used in this field.

Based on this framework, this project develops `MAXAID`, a lightweight `MATLAB-based` reimplementation. The code follows the original program logic and adopts a concise single-file structure, emphasizing readability and ease of use. It enables rapid implementation, testing, and comparison of different nonadiabatic models and algorithmic improvements. Functionally, `MAXAID` retains the electron–nuclear coupling formalism and incorporates improved surface hopping methods, including the `SDM` approach. Benefiting from `MATLAB`’s visualization and interactive capabilities, it is particularly suitable for method development, teaching, and prototyping, while supporting deployment across multiple platforms.

## Usage

Run `namd.m` in `MATLAB` or execute it via the command line using `matlab namd.m`.

## Citation

Original `PYXAID` references:

* Akimov A V, Prezhdo O V. The PYXAID program for non-adiabatic molecular dynamics in condensed matter systems. Journal of Chemical Theory and Computation, 2013, 9(11): 4959–4972.

* Akimov A V, Prezhdo O V. Advanced capabilities of the PYXAID program: integration schemes, decoherence effects, multiexcitonic states, and field–matter interaction. Journal of Chemical Theory and Computation, 2014, 10(2): 789–804.

This work:

To be added after the paper is officially published.