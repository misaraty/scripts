## **[中文版本](https://www.misaraty.com/2025-06-26_md%E9%94%AE%E9%95%BF%E9%94%AE%E8%A7%92%E5%88%86%E6%9E%90/)**

## **MD Simulation**

* **POSCAR**

Using [CaTiO₃ mp-5827](https://legacy.materialsproject.org/materials/mp-5827/) as an example, perform a high-precision SCF calculation to obtain a `POSCAR` file for `MD` input.

* **INCAR**

```shell
PREC = Accurate
SMASS = -3   #micro canonical ensemble
LREAL = A   #real space
IBRION = 0   #molecular dynamics
NBLOCK = 1   #update XDATCAR every x steps
TEBEG = 300   #start temperature
TEEND = 300   #final temperature
ISIF = 2   #2, ions change; 3, shape and ions change
ISYM = 0   #symmetry
NSW = 1000   #max ionic steps
POTIM = 1   #step width scaling
EDIFFG = -0.01   #ionic relaxation
EDIFF = 1e-5   #electronic SC-loop
LCHARG = .F.   #not save CHGCAR CHG
LWAVE = .F.   #not save WAVECAR
ISMEAR = 0   #gaussian smearing
SIGMA = 0.05   #the width of the smearing in eV
ALGO = Normal   #electronic minimisation algorithm
IVDW = 12   #DFT-D3 method with Becke-Jonson damping, van der Waals
```

> [!NOTE]
> If needed, use `SMASS = -1` for initial thermalization, then switch to `SMASS = -3` for temperature control.

### **Generating pfiles**

### **Method 1: Based on OUTCAR + CONTCAR**

Required: `OUTCAR` and `CONTCAR` from an `MD` run.

> [!NOTE]
> You may also use `OUTCAR + POSCAR` with slight modifications to the script.

Use the script `md_outcar_pfiles_v2.2.py` to generate pfiles.




















