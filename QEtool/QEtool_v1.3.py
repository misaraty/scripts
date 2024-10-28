import os
os.chdir(os.path.split(os.path.realpath(__file__))[0])
print('copyright by Zhaosheng Zhang (misaraty@163.com)\n' + 'last update: 2024-09-26\n')

'''
pseudopotentials:
  https://www.materialscloud.org/discover/sssp version 1.3.0
  https://dalcorso.github.io/pslibrary/ version 1.0.0
'''

upf_type = 'pslibrary'  # 'sssp_precision' | 'sssp_efficiency' | 'pslibrary'
calculation_type = 'dos_plot' # 'relax' | 'vc-relax' | 'scf' | 'nscf' | 'dos_plot'
poscar = 'cdte.vasp'
qe_bin = '/opt/ohpc/pub/apps/qe-7.3/bin'

elements = {
'H' : 1.00794,
'He' : 4.002602,
'Li' : 6.941,
'Be' : 9.012182,
'B' : 10.811,
'C' : 12.0107,
'N' : 14.0067,
'O' : 15.9994,
'F' : 18.9984032,
'Ne' : 20.1797,
'Na' : 22.98977,
'Mg' : 24.305,
'Al' : 26.981538,
'Si' : 28.0855,
'P' : 30.973761,
'S' : 32.065,
'Cl' : 35.453,
'Ar' : 39.948,
'K' : 39.0983,
'Ca' : 40.078,
'Sc' : 44.95591,
'Ti' : 47.867,
'V' : 50.9415,
'Cr' : 51.9961,
'Mn' : 54.938049,
'Fe' : 55.845,
'Co' : 58.9332,
'Ni' : 58.6934,
'Cu' : 63.546,
'Zn' : 65.39,
'Ga' : 69.723,
'Ge' : 72.64,
'As' : 74.9216,
'Se' : 78.96,
'Br' : 79.904,
'Kr' : 83.80,
'Rb' : 85.4678,
'Sr' : 87.62,
'Y' : 88.90585,
'Zr' : 91.224,
'Nb' : 92.90638,
'Mo' : 95.94,
'Tc' : 98,
'Ru' : 101.07,
'Rh' : 102.9055,
'Pd' : 106.42,
'Ag' : 107.8682,
'Cd' : 112.411,
'In' : 114.818,
'Sn' : 118.71,
'Sb' : 121.76,
'Te' : 127.6,
'I' : 126.90447,
'Xe' : 131.293,
'Cs' : 132.90545,
'Ba' : 137.327,
'La' : 138.9055,
'Ce' : 140.116,
'Pr' : 140.90765,
'Nd' : 144.24,
'Pm' : 145,
'Sm' : 150.36,
'Eu' : 151.964,
'Gd' : 157.25,
'Tb' : 158.92534,
'Dy' : 162.5,
'Ho' : 164.93032,
'Er' : 167.259,
'Tm' : 168.93421,
'Yb' : 173.04,
'Lu' : 174.967,
'Hf' : 178.49,
'Ta' : 180.9479,
'W' : 183.84,
'Re' : 186.207,
'Os' : 190.23,
'Ir' : 192.217,
'Pt' : 195.078,
'Au' : 196.96655,
'Hg' : 200.59,
'Tl' : 204.3833,
'Pb' : 207.2,
'Bi' : 208.98038,
'Po' : 209,
'At' : 210,
'Rn' : 222,
'Fr' : 223,
'Ra' : 226,
'Ac' : 227,
'Th' : 232.0381,
'Pa' : 231.03588,
'U' : 238.02891,
'Np' : 237,
'Pu' : 244,
'Am' : 243,
'Cm' : 247,
'Bk' : 247,
'Cf' : 251,
'Es' : 252,
'Fm' : 257,
'Md' : 258,
'No' : 259,
'Lr' : 262,
'Rf' : 261,
'Db' : 262,
'Sg' : 266.12,
'Bh' : 264.12,
'Hs' : 269.13,
'Mt' : 268.14,
'Ds' : 281.15,
'Rg' : 272.15,
'Cn' : 285,
}

upf_files_sssp_precision = """
Zr Zr_pbe_v1.uspp.F.UPF
B B_pbe_v1.01.uspp.F.UPF
In In.pbe-dn-rrkjus_psl.0.2.2.UPF
Na Na.paw.z_9.ld1.psl.v1.0.0-low.upf
Sr Sr_pbe_v1.uspp.F.UPF
Rh Rh_ONCV_PBE-1.0.oncvpsp.upf
V v_pbe_v1.4.uspp.F.UPF
Re Re_pbe_v1.2.uspp.F.UPF
Sn Sn_pbe_v1.uspp.F.UPF
Pt Pt.pbe-spfn-rrkjus_psl.1.0.0.UPF
Y Y_pbe_v1.uspp.F.UPF
Th Th.paw.z_12.ld1.uni-marburg.v0.upf
As As.nc.z_15.oncvpsp3.dojo.v4-std.upf
Cd Cd.paw.z_20.ld1.psl.v1.0.0-high.upf
Ni ni_pbe_v1.4.uspp.F.UPF
Fr Fr.paw.z_19.ld1.psl.v1.0.0-high.upf
Gd Gd.paw.z_18.atompaw.wentzcovitch.v1.2.upf
Mo Mo_ONCV_PBE-1.0.oncvpsp.upf
Ce Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf
Pm Pm.paw.z_15.atompaw.wentzcovitch.v1.2.upf
Pa Pa.paw.z_13.ld1.uni-marburg.v0.upf
Er Er.paw.z_22.atompaw.wentzcovitch.v1.2.upf
Np Np.paw.z_15.ld1.uni-marburg.v0.upf
Ir Ir.us.z_31.ld1.psl.v1.0.0-high.upf
Ac Ac.us.z_11.ld1.psl.v1.0.0-high.upf
Kr Kr.paw.z_18.ld1.psl.v1.0.0-high.upf
Cm Cm.paw.z_18.ld1.uni-marburg.v0.upf
Nd Nd.paw.z_14.atompaw.wentzcovitch.v1.2.upf
Br br_pbe_v1.4.uspp.F.UPF
Al Al.pbe-n-kjpaw_psl.1.0.0.UPF
F F.oncvpsp.upf
Mg mg_pbe_v1.4.uspp.F.UPF
Ag Ag_ONCV_PBE-1.0.oncvpsp.upf
Lr Lr.paw.z_25.ld1.uni-marburg.v0.upf
Tm Tm.paw.z_23.atompaw.wentzcovitch.v1.2.upf
Ge ge_pbe_v1.4.uspp.F.UPF
C C.pbe-n-kjpaw_psl.1.0.0.UPF
S s_pbe_v1.4.uspp.F.UPF
Dy Dy.paw.z_20.atompaw.wentzcovitch.v1.2.upf
Cf Cf.paw.z_20.ld1.uni-marburg.v0.upf
Pu Pu.paw.z_16.ld1.uni-marburg.v0.upf
Ar Ar.paw.z_8.ld1.psl.v1.0.0-high.upf
Tb Tb.paw.z_19.atompaw.wentzcovitch.v1.2.upf
Fm Fm.paw.z_22.ld1.uni-marburg.v0.upf
Es Es.paw.z_21.ld1.uni-marburg.v0.upf
Co Co_pbe_v1.2.uspp.F.UPF
Rn Rn.paw.z_18.ld1.psl.v1.0.0-high.upf
He He_ONCV_PBE-1.0.oncvpsp.upf
O O.pbe-n-kjpaw_psl.0.1.UPF
Bk Bk.paw.z_19.ld1.uni-marburg.v0.upf
No No.paw.z_24.ld1.uni-marburg.v0.upf
Ra Ra.paw.z_20.ld1.psl.v1.0.0-high.upf
Cr cr_pbe_v1.5.uspp.F.UPF
Hg Hg.us.z_12.uspp.gbrv.v1.upf
H H_ONCV_PBE-1.0.oncvpsp.upf
Eu Eu.paw.z_17.atompaw.wentzcovitch.v1.2.upf
Nb Nb.pbe-spn-kjpaw_psl.0.3.0.UPF
Pr Pr.paw.z_13.atompaw.wentzcovitch.v1.2.upf
Ru Ru_ONCV_PBE-1.0.oncvpsp.upf
Li li_pbe_v1.4.uspp.F.UPF
Tl Tl_pbe_v1.2.uspp.F.UPF
Te Te.us.z_6.ld1.psl.v1.0.0-low.upf
U U.paw.z_14.ld1.uni-marburg.v0.upf
N N.oncvpsp.upf
Rb Rb_ONCV_PBE-1.0.oncvpsp.upf
At At.us.z_17.ld1.psl.v1.0.0-high.upf
P P.pbe-n-rrkjus_psl.1.0.0.UPF
Ti ti_pbe_v1.4.uspp.F.UPF
Cu Cu.paw.z_11.ld1.psl.v1.0.0-low.upf
Be Be_ONCV_PBE-1.0.oncvpsp.upf
Au Au_ONCV_PBE-1.0.oncvpsp.upf
Zn Zn_pbe_v1.uspp.F.UPF
Po Po.pbe-dn-rrkjus_psl.1.0.0.UPF
Mn mn_pbe_v1.5.uspp.F.UPF
Sm Sm.paw.z_16.atompaw.wentzcovitch.v1.2.upf
Pb Pb.pbe-dn-kjpaw_psl.0.2.2.UPF
Ba Ba.nc.z_10.oncvpsp4.dojo.v4-sp.upf
La La.paw.z_11.atompaw.wentzcovitch.v1.2.upf
Ho Ho.paw.z_21.atompaw.wentzcovitch.v1.2.upf
Yb Yb.paw.z_24.atompaw.wentzcovitch.v1.2.upf
Xe Xe.paw.z_18.ld1.psl.v1.0.0-high.upf
Sc Sc.pbe-spn-kjpaw_psl.0.2.3.UPF
Cs Cs.nc.z_9.oncvpsp3.dojo.v4-str.upf
Ta Ta_pbe_v1.uspp.F.UPF
Cl Cl.pbe-n-rrkjus_psl.1.0.0.UPF
K K.pbe-spn-kjpaw_psl.1.0.0.UPF
Ne Ne.paw.z_8.ld1.psl.v1.0.0-high.upf
Pd Pd_ONCV_PBE-1.0.oncvpsp.upf
Md Md.paw.z_23.ld1.uni-marburg.v0.upf
Ga Ga.pbe-dn-kjpaw_psl.1.0.0.UPF
Ca Ca_pbe_v1.uspp.F.UPF
Sb sb_pbe_v1.4.uspp.F.UPF
I I.nc.z_17.oncvpsp4.sg15.v0.upf
Lu Lu.paw.z_25.atompaw.wentzcovitch.v1.2.upf
Am Am.paw.z_17.ld1.uni-marburg.v0.upf
Bi Bi_pbe_v1.uspp.F.UPF
W W_pbe_v1.2.uspp.F.UPF
Si Si.pbe-n-rrkjus_psl.1.0.0.UPF
Tc Tc_ONCV_PBE-1.0.oncvpsp.upf
Os Os_pbe_v1.2.uspp.F.UPF
Hf Hf-sp.oncvpsp.upf
Fe Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
Se Se_pbe_v1.uspp.F.UPF
""".strip()

upf_files_sssp_efficiency = """
Zr Zr_pbe_v1.uspp.F.UPF
Ar Ar_ONCV_PBE-1.1.oncvpsp.upf
Ir Ir_pbe_v1.2.uspp.F.UPF
In In.pbe-dn-rrkjus_psl.0.2.2.UPF
Cs Cs_pbe_v1.uspp.F.UPF
Sr Sr_pbe_v1.uspp.F.UPF
Hg Hg_ONCV_PBE-1.0.oncvpsp.upf
As As.pbe-n-rrkjus_psl.0.2.UPF
Rh Rh_ONCV_PBE-1.0.oncvpsp.upf
V v_pbe_v1.4.uspp.F.UPF
Re Re_pbe_v1.2.uspp.F.UPF
Sn Sn_pbe_v1.uspp.F.UPF
Y Y_pbe_v1.uspp.F.UPF
Th Th.paw.z_12.ld1.uni-marburg.v0.upf
Ni ni_pbe_v1.4.uspp.F.UPF
Fr Fr.paw.z_19.ld1.psl.v1.0.0-high.upf
Gd Gd.paw.z_18.atompaw.wentzcovitch.v1.2.upf
Mo Mo_ONCV_PBE-1.0.oncvpsp.upf
Ce Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf
Pm Pm.paw.z_15.atompaw.wentzcovitch.v1.2.upf
Pa Pa.paw.z_13.ld1.uni-marburg.v0.upf
Er Er.paw.z_22.atompaw.wentzcovitch.v1.2.upf
Np Np.paw.z_15.ld1.uni-marburg.v0.upf
Te Te_pbe_v1.uspp.F.UPF
Ac Ac.us.z_11.ld1.psl.v1.0.0-high.upf
Cm Cm.paw.z_18.ld1.uni-marburg.v0.upf
Nd Nd.paw.z_14.atompaw.wentzcovitch.v1.2.upf
Br br_pbe_v1.4.uspp.F.UPF
Al Al.pbe-n-kjpaw_psl.1.0.0.UPF
Ag Ag_ONCV_PBE-1.0.oncvpsp.upf
Lr Lr.paw.z_25.ld1.uni-marburg.v0.upf
Tm Tm.paw.z_23.atompaw.wentzcovitch.v1.2.upf
B b_pbe_v1.4.uspp.F.UPF
Sc Sc_ONCV_PBE-1.0.oncvpsp.upf
Be be_pbe_v1.4.uspp.F.UPF
Ge ge_pbe_v1.4.uspp.F.UPF
C C.pbe-n-kjpaw_psl.1.0.0.UPF
S s_pbe_v1.4.uspp.F.UPF
Dy Dy.paw.z_20.atompaw.wentzcovitch.v1.2.upf
Cf Cf.paw.z_20.ld1.uni-marburg.v0.upf
Pu Pu.paw.z_16.ld1.uni-marburg.v0.upf
Tb Tb.paw.z_19.atompaw.wentzcovitch.v1.2.upf
Kr Kr_ONCV_PBE-1.0.oncvpsp.upf
Rn Rn.pbe-dn-kjpaw_psl.1.0.0.UPF
Fm Fm.paw.z_22.ld1.uni-marburg.v0.upf
Es Es.paw.z_21.ld1.uni-marburg.v0.upf
Co Co_pbe_v1.2.uspp.F.UPF
He He_ONCV_PBE-1.0.oncvpsp.upf
O O.pbe-n-kjpaw_psl.0.1.UPF
Bk Bk.paw.z_19.ld1.uni-marburg.v0.upf
No No.paw.z_24.ld1.uni-marburg.v0.upf
Ra Ra.paw.z_20.ld1.psl.v1.0.0-high.upf
Cr cr_pbe_v1.5.uspp.F.UPF
I I.pbe-n-kjpaw_psl.0.2.UPF
Eu Eu.paw.z_17.atompaw.wentzcovitch.v1.2.upf
Nb Nb.pbe-spn-kjpaw_psl.0.3.0.UPF
Mg Mg.pbe-n-kjpaw_psl.0.3.0.UPF
Pr Pr.paw.z_13.atompaw.wentzcovitch.v1.2.upf
Ru Ru_ONCV_PBE-1.0.oncvpsp.upf
Li li_pbe_v1.4.uspp.F.UPF
Tl Tl_pbe_v1.2.uspp.F.UPF
U U.paw.z_14.ld1.uni-marburg.v0.upf
Rb Rb_ONCV_PBE-1.0.oncvpsp.upf
At At.us.z_17.ld1.psl.v1.0.0-high.upf
P P.pbe-n-rrkjus_psl.1.0.0.UPF
Ti ti_pbe_v1.4.uspp.F.UPF
Xe Xe_ONCV_PBE-1.1.oncvpsp.upf
Cu Cu.paw.z_11.ld1.psl.v1.0.0-low.upf
Pt pt_pbe_v1.4.uspp.F.UPF
Na na_pbe_v1.5.uspp.F.UPF
Cd Cd.pbe-dn-rrkjus_psl.0.3.1.UPF
Au Au_ONCV_PBE-1.0.oncvpsp.upf
Zn Zn_pbe_v1.uspp.F.UPF
Po Po.pbe-dn-rrkjus_psl.1.0.0.UPF
Ne Ne_ONCV_PBE-1.0.oncvpsp.upf
F f_pbe_v1.4.uspp.F.UPF
N N.pbe-n-radius_5.UPF
Mn mn_pbe_v1.5.uspp.F.UPF
Sm Sm.paw.z_16.atompaw.wentzcovitch.v1.2.upf
Pb Pb.pbe-dn-kjpaw_psl.0.2.2.UPF
La La.paw.z_11.atompaw.wentzcovitch.v1.2.upf
Ho Ho.paw.z_21.atompaw.wentzcovitch.v1.2.upf
Yb Yb.paw.z_24.atompaw.wentzcovitch.v1.2.upf
Ba Ba.pbe-spn-kjpaw_psl.1.0.0.UPF
Cl cl_pbe_v1.4.uspp.F.UPF
Ta Ta_pbe_v1.uspp.F.UPF
K K.pbe-spn-kjpaw_psl.1.0.0.UPF
Pd Pd_ONCV_PBE-1.0.oncvpsp.upf
Md Md.paw.z_23.ld1.uni-marburg.v0.upf
Ca Ca_pbe_v1.uspp.F.UPF
Ga Ga.pbe-dn-kjpaw_psl.1.0.0.UPF
H H.pbe-rrkjus_psl.1.0.0.UPF
Sb sb_pbe_v1.4.uspp.F.UPF
Lu Lu.paw.z_25.atompaw.wentzcovitch.v1.2.upf
Am Am.paw.z_17.ld1.uni-marburg.v0.upf
Bi Bi_pbe_v1.uspp.F.UPF
W W_pbe_v1.2.uspp.F.UPF
Si Si.pbe-n-rrkjus_psl.1.0.0.UPF
Tc Tc_ONCV_PBE-1.0.oncvpsp.upf
Os Os_pbe_v1.2.uspp.F.UPF
Hf Hf-sp.oncvpsp.upf
Fe Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
Se Se_pbe_v1.uspp.F.UPF
""".strip()

upf_files_pslibrary = """
Ac Ac.pbe-spfn-kjpaw_psl.1.0.0.UPF
Ag Ag.pbe-n-kjpaw_psl.1.0.0.UPF
Al Al.pbe-n-kjpaw_psl.1.0.0.UPF
Am Am.pbe-spfn-kjpaw_psl.1.0.0.UPF
Ar Ar.pbe-n-kjpaw_psl.1.0.0.UPF
As As.pbe-n-kjpaw_psl.1.0.0.UPF
At At.pbe-dn-kjpaw_psl.1.0.0.UPF
Au Au.pbe-n-kjpaw_psl.1.0.0.UPF
B B.pbe-n-kjpaw_psl.1.0.0.UPF
Ba Ba.pbe-spn-kjpaw_psl.1.0.0.UPF
Be Be.pbe-n-kjpaw_psl.1.0.0.UPF
Bi Bi.pbe-dn-kjpaw_psl.1.0.0.UPF
Br Br.pbe-n-kjpaw_psl.1.0.0.UPF
C C.pbe-n-kjpaw_psl.1.0.0.UPF
Ca Ca.pbe-spn-kjpaw_psl.1.0.0.UPF
Cd Cd.pbe-n-kjpaw_psl.1.0.0.UPF
Ce Ce.pbe-spdn-kjpaw_psl.1.0.0.UPF
Cl Cl.pbe-n-kjpaw_psl.1.0.0.UPF
Co Co.pbe-n-kjpaw_psl.1.0.0.UPF
Cr Cr.pbe-spn-kjpaw_psl.1.0.0.UPF
Cs Cs.pbe-spn-kjpaw_psl.1.0.0.UPF
Cu Cu.pbe-dn-kjpaw_psl.1.0.0.UPF
Dy Dy.pbe-spdn-kjpaw_psl.1.0.0.UPF
Er Er.pbe-spdn-kjpaw_psl.1.0.0.UPF
F F.pbe-n-kjpaw_psl.1.0.0.UPF
Eu Eu.pbe-spn-kjpaw_psl.1.0.0.UPF
Fe Fe.pbe-n-kjpaw_psl.1.0.0.UPF
Fr Fr.pbe-spdn-kjpaw_psl.1.0.0.UPF
Ga Ga.pbe-dn-kjpaw_psl.1.0.0.UPF
Gd Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF
Ge Ge.pbe-n-kjpaw_psl.1.0.0.UPF
H H.pbe-kjpaw_psl.1.0.0.UPF
He He.pbe-kjpaw_psl.1.0.0.UPF
Hf Hf.pbe-spn-kjpaw_psl.1.0.0.UPF
Hg Hg.pbe-n-kjpaw_psl.1.0.0.UPF
Ho Ho.pbe-spdn-kjpaw_psl.1.0.0.UPF
I. I.pbe-n-kjpaw_psl.1.0.0.UPF
In In.pbe-dn-kjpaw_psl.1.0.0.UPF
Ir Ir.pbe-n-kjpaw_psl.1.0.0.UPF
K K.pbe-spn-kjpaw_psl.1.0.0.UPF
Kr Kr.pbe-dn-kjpaw_psl.1.0.0.UPF
La La.pbe-spfn-kjpaw_psl.1.0.0.UPF
Li Li.pbe-s-kjpaw_psl.1.0.0.UPF
Lu Lu.pbe-spdn-kjpaw_psl.1.0.0.UPF
Mg Mg.pbe-spn-kjpaw_psl.1.0.0.UPF
Mn Mn.pbe-spn-kjpaw_psl.1.0.0.UPF
Mo Mo.pbe-spn-kjpaw_psl.1.0.0.UPF
N. N.pbe-n-kjpaw_psl.1.0.0.UPF
Na Na.pbe-spn-kjpaw_psl.1.0.0.UPF
Nb Nb.pbe-spn-kjpaw_psl.1.0.0.UPF
Ne Ne.pbe-n-kjpaw_psl.1.0.0.UPF
Ni Ni.pbe-n-kjpaw_psl.1.0.0.UPF
Np Np.pbe-spfn-kjpaw_psl.1.0.0.UPF
O O.pbe-n-kjpaw_psl.1.0.0.UPF
Os Os.pbe-spn-kjpaw_psl.1.0.0.UPF
P P.pbe-n-kjpaw_psl.1.0.0.UPF
Pa Pa.pbe-spfn-kjpaw_psl.1.0.0.UPF
Pb Pb.pbe-dn-kjpaw_psl.1.0.0.UPF
Pd Pd.pbe-n-kjpaw_psl.1.0.0.UPF
Pm Pm.pbe-spdn-kjpaw_psl.1.0.0.UPF
Po Po.pbe-dn-kjpaw_psl.1.0.0.UPF
Pr Pr.pbe-spdn-kjpaw_psl.1.0.0.UPF
Pt Pt.pbe-n-kjpaw_psl.1.0.0.UPF
Pu Pu.pbe-spfn-kjpaw_psl.1.0.0.UPF
Ra Ra.pbe-spdn-kjpaw_psl.1.0.0.UPF
Rb Rb.pbe-spn-kjpaw_psl.1.0.0.UPF
Re Re.pbe-spn-kjpaw_psl.1.0.0.UPF
Rh Rh.pbe-spn-kjpaw_psl.1.0.0.UPF
Rn Rn.pbe-dn-kjpaw_psl.1.0.0.UPF
Ru Ru.pbe-spn-kjpaw_psl.1.0.0.UPF
S S.pbe-n-kjpaw_psl.1.0.0.UPF
Sb Sb.pbe-n-kjpaw_psl.1.0.0.UPF
Sc Sc.pbe-spn-kjpaw_psl.1.0.0.UPF
Se Se.pbe-n-kjpaw_psl.1.0.0.UPF
Si Si.pbe-n-kjpaw_psl.1.0.0.UPF
Sm Sm.pbe-spdn-kjpaw_psl.1.0.0.UPF
Sn Sn.pbe-dn-kjpaw_psl.1.0.0.UPF
Sr Sr.pbe-spn-kjpaw_psl.1.0.0.UPF
Ta Ta.pbe-spn-kjpaw_psl.1.0.0.UPF
Tb Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF
Tc Tc.pbe-spn-kjpaw_psl.1.0.0.UPF
Te Te.pbe-n-kjpaw_psl.1.0.0.UPF
Th Th.pbe-spfn-kjpaw_psl.1.0.0.UPF
Ti Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
Tl Tl.pbe-dn-kjpaw_psl.1.0.0.UPF
Tm Tm.pbe-spdn-kjpaw_psl.1.0.0.UPF
U U.pbe-spfn-kjpaw_psl.1.0.0.UPF
V V.pbe-spn-kjpaw_psl.1.0.0.UPF
W W.pbe-spn-kjpaw_psl.1.0.0.UPF
Xe Xe.pbe-dn-kjpaw_psl.1.0.0.UPF
Y Y.pbe-spn-kjpaw_psl.1.0.0.UPF
Yb Yb.pbe-spn-kjpaw_psl.1.0.0.UPF
Zn Zn.pbe-dn-kjpaw_psl.1.0.0.UPF
Zr Zr.pbe-spn-kjpaw_psl.1.0.0.UPF
""".strip()

if upf_type == 'sssp_precision':
    upf_files = {line.split()[0]: line.split()[1] for line in upf_files_sssp_precision.strip().split('\n')}
    pseudo_dir = '/opt/ohpc/pub/apps/sssp/precision/'
elif upf_type == 'sssp_efficiency':
    upf_files = {line.split()[0]: line.split()[1] for line in upf_files_sssp_efficiency.strip().split('\n')}
    pseudo_dir = '/opt/ohpc/pub/apps/sssp/efficiency/'
elif upf_type == 'pslibrary':
    upf_files = {line.split()[0]: line.split()[1] for line in upf_files_pslibrary.strip().split('\n')}
    pseudo_dir = '/opt/ohpc/pub/apps/pslibrary.1.0.0/pbe/PSEUDOPOTENTIALS/'

with open(f'{poscar}', 'r') as file:
    poscar_lines = file.readlines()

lattice_vectors = poscar_lines[2:5]
atom_types = poscar_lines[5].split()
num_atoms_per_type = list(map(int, poscar_lines[6].split()))
total_atoms = sum(num_atoms_per_type)

cell_parameters = 'CELL_PARAMETERS angstrom\n' + ''.join(lattice_vectors)
atomic_species = 'ATOMIC_SPECIES\n'
atomic_positions = 'ATOMIC_POSITIONS crystal\n'
start_line = 8

for atom_type, num_atoms in zip(atom_types, num_atoms_per_type):
    atomic_species += f"{atom_type}     {elements[atom_type]}     {upf_files.get(atom_type, 'NOT_FOUND.UPF')}\n"
    for i in range(start_line, start_line + num_atoms):
        atomic_positions += f"{atom_type}     {poscar_lines[i].strip()}\n"
    start_line += num_atoms

atomic_species = atomic_species.strip()
atomic_positions = atomic_positions.strip()

if calculation_type == 'relax':
    qe_input = f"""&CONTROL
  calculation = '{calculation_type}'
  outdir = './out/'
  prefix = 'qe'
  pseudo_dir = '{pseudo_dir}'
  nstep = 1000
/
&SYSTEM
  ecutwfc = 60
  ecutrho = 600
  ibrav = 0
  nat = {total_atoms}
  ntyp = {len(atom_types)}
  vdw_corr = 'DFT-D3'
/
&ELECTRONS
  conv_thr = 1.D-7
  electron_maxstep = 1000
  mixing_beta = 0.7
/
&IONS
  ion_dynamics = 'bfgs'
/
&CELL
  cell_dynamics = 'bfgs'
/
{atomic_species}
{atomic_positions}
K_POINTS automatic
4 4 4 0 0 0
{cell_parameters}
"""

elif calculation_type == 'vc-relax':
    qe_input = f"""&CONTROL
  calculation = '{calculation_type}'
  outdir = './out/'
  prefix = 'qe'
  pseudo_dir = '{pseudo_dir}'
  nstep = 1000
/
&SYSTEM
  ecutwfc = 60
  ecutrho = 600
  ibrav = 0
  nat = {total_atoms}
  ntyp = {len(atom_types)}
  vdw_corr = 'DFT-D3'
/
&ELECTRONS
  conv_thr = 1.D-7
  electron_maxstep = 1000
  mixing_beta = 0.7
/
&IONS
  ion_dynamics = 'bfgs'
/
&CELL
  cell_dynamics = 'bfgs'
  cell_dofree = 'all'
/
{atomic_species}
{atomic_positions}
K_POINTS automatic
4 4 4 0 0 0
{cell_parameters}
"""

elif calculation_type == 'scf':

    if os.path.exists('vc-relax.out'):
        with open('vc-relax.out', 'r') as file:
            lines = file.readlines()
        start = len(lines) - lines[::-1].index("ATOMIC_POSITIONS (crystal)\n")
        end = lines.index("End final coordinates\n")
        atomic_positions = 'ATOMIC_POSITIONS crystal\n' + ''.join(lines[start:end]).strip()
        cell_parameters = 'CELL_PARAMETERS angstrom\n' + ''.join(lines[start-5:start-1])

    elif os.path.exists('relax.out'):
        with open('relax.out', 'r') as file:
            lines = file.readlines()
        start = lines.index("ATOMIC_POSITIONS (crystal)\n") + 1
        end = lines.index("End final coordinates\n")
        atomic_positions = 'ATOMIC_POSITIONS crystal\n' + ''.join(lines[start:end]).strip()

    qe_input = f"""&CONTROL
  calculation = '{calculation_type}'
  outdir = './out/'
  prefix = 'qe'
  pseudo_dir = '{pseudo_dir}'
/
&SYSTEM
  ecutwfc = 60
  ecutrho = 600
  ibrav = 0
  nat = {total_atoms}
  ntyp = {len(atom_types)}
  occupations = 'smearing'
  smearing = 'gaussian'
  degauss = 0.001
/
&ELECTRONS
  conv_thr = 1.D-7
  electron_maxstep = 1000
  mixing_beta = 0.7
/
{atomic_species}
{atomic_positions}
K_POINTS automatic
8 8 8 0 0 0
{cell_parameters}
"""

elif calculation_type == 'nscf':

    if os.path.exists('vc-relax.out'):
        with open('vc-relax.out', 'r') as file:
            lines = file.readlines()
        start = len(lines) - lines[::-1].index("ATOMIC_POSITIONS (crystal)\n")
        end = lines.index("End final coordinates\n")
        atomic_positions = 'ATOMIC_POSITIONS crystal\n' + ''.join(lines[start:end]).strip()
        cell_parameters = 'CELL_PARAMETERS angstrom\n' + ''.join(lines[start-5:start-1])

    elif os.path.exists('relax.out'):
        with open('relax.out', 'r') as file:
            lines = file.readlines()
        start = lines.index("ATOMIC_POSITIONS (crystal)\n") + 1
        end = lines.index("End final coordinates\n")
        atomic_positions = 'ATOMIC_POSITIONS crystal\n' + ''.join(lines[start:end]).strip()

    qe_input = f"""&CONTROL
  calculation = '{calculation_type}'
  outdir = './out/'
  prefix = 'qe'
  pseudo_dir = '{pseudo_dir}'
/
&SYSTEM
  ecutwfc = 60
  ecutrho = 600
  ibrav = 0
  nat = {total_atoms}
  ntyp = {len(atom_types)}
  occupations = 'tetrahedra'
  degauss = 0.001
/
&ELECTRONS
  conv_thr = 1.D-7
  electron_maxstep = 1000
  mixing_beta = 0.7
/
{atomic_species}
{atomic_positions}
K_POINTS automatic
8 8 8 0 0 0
{cell_parameters}
"""

elif calculation_type == 'dos_plot':
    import subprocess
    import time
    import numpy as np
    import matplotlib.pyplot as plt

    with open('pdos.in', 'w') as f:
        f.write("&PROJWFC\n")
        f.write("    prefix = 'qe'\n")
        f.write("    outdir = './out/'\n")
        f.write("    Emin=-20.0, Emax=20.0, DeltaE=0.1,\n")
        f.write("    ngauss=0, degauss=0.01\n")
        f.write("/\n")

    subprocess.run(f"{qe_bin}/projwfc.x < pdos.in > pdos.out", shell=True)

    for element in atom_types:
        subprocess.run(f"{qe_bin}/sumpdos.x *\\({element}\\)* > {element}.dat", shell=True)

    with open('nscf.out', 'r') as f:
        for line in f:
            if "the Fermi energy is" in line:
                fermi_energy = float(line.split()[-2])

    dos_total = []
    pdos_dict = {element: [] for element in atom_types}

    with open('qe.pdos_tot', 'r') as f:
        for line in f:
            if not line.startswith('#'):
                columns = list(map(float, line.split()))
                dos_total.append((columns[0] - fermi_energy, columns[-1]))

    for element in atom_types:
        with open(f"{element}.dat", 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    columns = list(map(float, line.split()))
                    pdos_dict[element].append((columns[0] - fermi_energy, columns[-1]))

    dos_total = np.array(dos_total)
    tab_colors = plt.get_cmap('tab10')
    num_elements = len(pdos_dict)

    total_color = tab_colors(0)
    plt.plot(dos_total[:, 0], dos_total[:, 1], label='Total', color=total_color)
    plt.fill_between(dos_total[:, 0], dos_total[:, 1], color=total_color, alpha=0.3)

    for i, (element, pdos) in enumerate(pdos_dict.items()):
        pdos = np.array(pdos)
        color = tab_colors(i + 1)
        plt.plot(pdos[:, 0], pdos[:, 1], label=f'{element}', color=color)
        plt.fill_between(pdos[:, 0], pdos[:, 1], color=color, alpha=0.3)

    plt.xlim(-3, 3)
    plt.ylim(0, 15)
    plt.xlabel('Energy (eV)', fontsize=16)
    plt.ylabel('PDOS (1/eV)', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.axvline(0, color='grey', lw=1, ls='--')
    plt.legend(fontsize=16, loc='best')
    plt.tight_layout()
    plt.savefig('dos.jpg', dpi=300)
    plt.close()

if calculation_type == 'relax' or calculation_type == 'vc-relax' or calculation_type == 'scf' or calculation_type == 'nscf':
    with open(f'{calculation_type}.in', 'w') as output_file:
        output_file.write(qe_input.strip() + '\n')

    print(f"Quantum Espresso input file '{calculation_type}.in' has been created.")

