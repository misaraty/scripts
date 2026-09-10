## [中文版本](https://www.misaraty.com/2026-05-02_maxaid/)

## MAXAID

`PYXAID`, developed by `Oleg V. Prezhdo` and `Alexey V. Akimov`, is a well-established nonadiabatic molecular dynamics software widely used in excited-state simulations of condensed matter systems. `Libra` further extends this framework. In addition, other software packages, such as `Hefei-NAMD`, `NEXMD`, `SHARC`, `Newton-X`, and `SPADE`, are also widely used in this field.

Based on `PYXAID`, this project develops `MAXAID`, a lightweight `MATLAB-based` reimplementation. The code follows the original program logic and adopts a concise single-file structure, emphasizing readability and ease of use. It enables rapid implementation, testing, and comparison of different nonadiabatic models and algorithmic improvements. Functionally, `MAXAID` retains the electron–nuclear coupling formalism and incorporates improved surface hopping methods, including the `SDM` approach. Benefiting from `MATLAB`’s visualization and interactive capabilities, it is particularly suitable for method development, teaching, and prototyping, while supporting deployment across multiple platforms.

## Usage

Run `namd.m` in `MATLAB` or execute it via the command line using `matlab namd.m`.

## Citation

Original `PYXAID` references:

* [Akimov, A. V.; Prezhdo, O. V*. The Pyxaid Program for Non-Adiabatic Molecular Dynamics in Condensed Matter Systems. J. Chem. Theory Comput. 2013, 9, 4959–4972](https://pubs.acs.org/jctcce/article-abstract/9/11/4959/843563/The-PYXAID-Program-for-Non-Adiabatic-Molecular?redirectedFrom=fulltext)

* [Akimov, A. V.; Prezhdo, O. V*. Advanced Capabilities of the Pyxaid Program: Integration Schemes, Decoherence Effects, Multiexcitonic States, and Field-Matter Interaction. J. Chem. Theory Comput. 2014, 10, 789–804](https://pubs.acs.org/jctcce/article/10/2/789/794232/Advanced-Capabilities-of-the-PYXAID-Program)

This work:

* [Zhang, Z.*; Liu, Y.; Liu, J. Phosphonic Acid Molecular Regulation of Frenkel Defects for Suppressing Nonradiative Recombination in FAPbI3 Perovskites. J. Phys. Chem. Lett. 2026, 17, 9756–9765](https://pubs.acs.org/jpclcd/article-abstract/17/33/9756/5250682/Phosphonic-Acid-Molecular-Regulation-of-Frenkel?redirectedFrom=fulltext)

* [Zhang, Z.*; Liu, J.; Liu, Y. Molecular Passivation of Iodine Vacancies Suppresses Nonradiative Recombination in FAPbI3 Perovskites. J. Mater. Chem. A 2026](https://pubs.rsc.org/ta/article-abstract/doi/10.1039/d6ta06718b/1300128/Molecular-passivation-of-iodine-vacancies?redirectedFrom=fulltext)