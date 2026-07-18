from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
RT = SCRIPT_DIR.parent
OUTDIR = SCRIPT_DIR / "out"
OUTDIR.mkdir(exist_ok=True)

namdtime = 10000
num_sh_traj = 1000
iconds_i = 1000
Ham_size = 4001

active_space = np.array([1, 2], dtype=int)
states_def = [
    ("GS", [1, -1], 0.00),
    ("S1", [1, -2], 0.00),
]

decoherence_method = 2
nucl_dt = 1.0
elec_dt = 1.0
integrator = 0
sh_algo = 0
boltz_flag = 1
Temp = 300.0
alp_bet = 0

hbar = 0.658218
kb = 8.617e-5
Ry_to_eV = 13.60569253
Ha_to_eV = 27.211396

active_map = {int(v): i for i, v in enumerate(active_space)}


@dataclass
class ElectronicStructure:
    num_states: int
    curr_state: int = 0
    Ccurr: np.ndarray = field(init=False)
    Cprev: np.ndarray = field(init=False)
    Cnext: np.ndarray = field(init=False)
    Hcurr: np.ndarray = field(init=False)
    Hprev: np.ndarray = field(init=False)
    Hnext: np.ndarray = field(init=False)
    dHdt: np.ndarray = field(init=False)
    g: np.ndarray = field(init=False)
    tau_m: np.ndarray = field(init=False)
    t_m: np.ndarray = field(init=False)
    A: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        n = self.num_states
        self.Ccurr = np.zeros(n, dtype=np.complex128)
        self.Cprev = np.zeros(n, dtype=np.complex128)
        self.Cnext = np.zeros(n, dtype=np.complex128)
        self.Hcurr = np.zeros((n, n), dtype=np.complex128)
        self.Hprev = np.zeros((n, n), dtype=np.complex128)
        self.Hnext = np.zeros((n, n), dtype=np.complex128)
        self.dHdt = np.zeros((n, n), dtype=np.complex128)
        self.g = np.eye(n, dtype=float)
        self.tau_m = np.zeros(n, dtype=float)
        self.t_m = np.zeros(n, dtype=float)
        self.A = np.zeros((n, n), dtype=np.complex128)


@dataclass
class MEState:
    name: str
    actual_state: np.ndarray
    shift: float = 0.0


def extract_2d(arr: np.ndarray, templ: np.ndarray, shift: int) -> np.ndarray:
    """MATLAB: idx = templ + shift + 1; out = in(idx, idx)."""
    idx = templ + shift
    return arr[np.ix_(idx, idx)]


def make_me_state(cs: Tuple[str, List[int], float]) -> MEState:
    name = cs[0]
    actual_state = np.array(cs[1], dtype=int)
    shift = float(cs[2]) if len(cs) >= 3 else 0.0
    return MEState(name=name, actual_state=actual_state, shift=shift)


def delta(
    A: np.ndarray, B: np.ndarray, a: int = 0, b: int = 0
) -> Tuple[bool, np.ndarray, np.ndarray, int, int]:
    C_ = np.unique(np.concatenate([A, B]))
    nexc = 0
    for val in C_:
        n_in_a = np.sum(A == val)
        n_in_b = np.sum(B == val)
        d = int(n_in_a - n_in_b)
        if d > 0:
            nexc += d
        if d == 1:
            a = int(val)
        if d == -1:
            b = int(val)
    return nexc == 1, A, B, a, b


def ext2int(external: int) -> int:
    idx = active_map[abs(int(external))]
    f = 1 if external < 0 else 0
    return 2 * abs(idx) + f


def regression(
    X: np.ndarray, Y: np.ndarray, a: float = 0.0, b: float = 0.0
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    sx = np.sum(X)
    if abs(sx) < np.finfo(float).eps:
        b = 0.0
    else:
        b = float(np.sum(Y) / sx)
    return X, Y, a, b


def set_state(es: ElectronicStructure, indx: int) -> ElectronicStructure:
    es.Ccurr[:] = 0.0
    es.Ccurr[indx] = 1.0 + 0.0j
    es.curr_state = int(indx)
    return es


def insertion(
    this: ElectronicStructure, es: ElectronicStructure
) -> ElectronicStructure:
    this.num_states = es.num_states
    this.curr_state = es.curr_state
    this.Cprev = es.Cprev.copy()
    this.Ccurr = es.Ccurr.copy()
    this.Cnext = es.Cnext.copy()
    this.tau_m = es.tau_m.copy()
    this.t_m = es.t_m.copy()
    this.A = es.A.copy()
    return this


def init_hop_prob1(es: ElectronicStructure) -> ElectronicStructure:
    es.g = np.eye(es.num_states, dtype=float)
    return es


def efield(t: float, E: np.ndarray, Eex: float) -> Tuple[np.ndarray, float]:
    return np.zeros_like(E), 0.0


def rot1(es: ElectronicStructure, phi: float, i: int, j: int) -> ElectronicStructure:
    c = math.cos(phi)
    s = math.sin(phi)
    ci = es.Ccurr[i]
    cj = es.Ccurr[j]
    es.Ccurr[i] = c * ci + s * cj
    es.Ccurr[j] = -s * ci + c * cj
    return es


def rot2(es: ElectronicStructure, phi: float, i: int, j: int) -> ElectronicStructure:
    cs = math.cos(phi)
    isi = 1j * math.sin(phi)
    ci = es.Ccurr[i]
    cj = es.Ccurr[j]
    es.Ccurr[i] = cs * ci + isi * cj
    es.Ccurr[j] = isi * ci + cs * cj
    return es


def rot(
    es: ElectronicStructure, Hij: complex, dt: float, i: int, j: int
) -> ElectronicStructure:
    phi1 = 0.5 * dt * Hij.imag / hbar
    phi2 = -dt * Hij.real / hbar
    es = rot1(es, phi1, i, j)
    es = rot2(es, phi2, i, j)
    es = rot1(es, phi1, i, j)
    return es


def phase(
    es: ElectronicStructure, Hii: complex, dt: float, i: int
) -> ElectronicStructure:
    phi = -dt * Hii.real / hbar
    es.Ccurr[i] = np.exp(1j * phi) * es.Ccurr[i]
    return es


def propagate_coefficients(
    es: ElectronicStructure, dt: float, Ef: np.ndarray
) -> ElectronicStructure:
    n = es.num_states
    for i in range(n):
        for j in range(i + 1, n):
            es = rot(es, es.Hcurr[i, j], 0.5 * dt, i, j)
    for i in range(n):
        es = phase(es, es.Hcurr[i, i], dt, i)
    for i in range(n - 1, -1, -1):
        for j in range(n - 1, i, -1):
            es = rot(es, es.Hcurr[i, j], 0.5 * dt, i, j)
    return es


def update_populations(es: ElectronicStructure) -> ElectronicStructure:
    es.A = np.outer(np.conj(es.Ccurr), es.Ccurr)
    return es


def update_hop_prob_fssh(
    es: ElectronicStructure,
    dt: float,
    boltz_flag_value: int,
    temp: float,
    Ef: np.ndarray,
    Eex: float,
    rates: np.ndarray,
) -> ElectronicStructure:
    es = update_populations(es)
    Heff = es.Hcurr
    n = es.num_states
    for i in range(n):
        a_ii = float(np.real(es.A[i, i]))
        if a_ii < 1e-12:
            a_ii = 1e-12
        g_row_sum = 0.0
        for j in range(n):
            if j != i:
                g_ij = (2.0 * dt / (a_ii * hbar)) * np.imag(es.A[i, j] * Heff[i, j])
                if g_ij < 0.0:
                    g_ij = 0.0
                E_i = float(np.real(Heff[i, i]))
                E_j = float(np.real(Heff[j, j]))
                dE = E_j - E_i
                bf = 1.0
                if dE > Eex:
                    bf = math.exp(-((dE - Eex) / (kb * temp)))
                g_ij *= bf
                g_row_sum += g_ij
                es.g[i, j] = g_ij
        es.g[i, i] = 1.0 - g_row_sum
    return es


def propagate_electronic(
    es_list: List[ElectronicStructure], i: int, rates: np.ndarray
) -> Tuple[List[ElectronicStructure], np.ndarray]:
    nel = int(round(nucl_dt / elec_dt))
    Eex = 0.0
    Ef = np.zeros(3, dtype=float)
    if integrator == 0:
        for j in range(nel):
            tim = i * nucl_dt + j * elec_dt
            Ef, Eex = efield(tim, Ef, Eex)
            es_list[i] = propagate_coefficients(es_list[i], elec_dt, Ef)
            if sh_algo == 0:
                es_list[i] = update_hop_prob_fssh(
                    es_list[i], elec_dt, boltz_flag, Temp, Ef, Eex, rates
                )
    return es_list, rates


def hop(sh_prob: np.ndarray, hopstate: int, numstates: int) -> Tuple[np.ndarray, int]:
    in_state = int(hopstate)
    ksi = np.random.rand()
    probs = sh_prob[in_state, :].copy()
    nrm = np.sum(probs)
    if nrm <= 0:
        return sh_prob, hopstate
    probs /= nrm
    cprobs = np.cumsum(probs)
    candidates = np.where((0 < cprobs) & (ksi <= cprobs))[0]
    hstate = int(candidates[0]) if candidates.size else in_state
    return sh_prob, hstate


def update_decoherence_times(
    es: ElectronicStructure, rates: np.ndarray
) -> ElectronicStructure:
    es = update_populations(es)
    pops = np.real(np.diag(es.A))
    es.tau_m = np.real(rates) @ pops
    return es


def decohere(es: ElectronicStructure, i: int) -> ElectronicStructure:
    es.Ccurr[:] = 0.0
    es.Ccurr[i] = 1.0
    es.curr_state = int(i)
    return es


def project_out(es: ElectronicStructure, i: int) -> ElectronicStructure:
    es = update_populations(es)
    es.Ccurr[i] = 0.0
    nrm2 = np.sum(np.real(np.diag(es.A))) - float(np.real(es.A[i, i]))
    if nrm2 <= 0:
        es.Ccurr[:] = 0.0
        es.Ccurr[es.curr_state] = 1.0
    else:
        es.Ccurr /= math.sqrt(nrm2)
    es = update_populations(es)
    return es


def dish_decoherence(
    es: ElectronicStructure,
    dt: float,
    boltz_flag_value: int,
    temp: float,
    rates: np.ndarray,
) -> ElectronicStructure:
    es = update_decoherence_times(es, rates)
    for i in range(es.num_states):
        if es.tau_m[i] <= 0:
            continue
        rnd_i = 1.0 / es.tau_m[i]
        if es.t_m[i] >= rnd_i:
            zeta = np.random.rand()
            P = float(np.real(es.A[i, i]))
            dE = float(np.real(es.Hcurr[i, i] - es.Hcurr[es.curr_state, es.curr_state]))
            if dE > 0:
                P *= math.exp(-(dE / (kb * temp)))
            if zeta < P:
                es = decohere(es, i)
                break
            es = project_out(es, i)
            es.t_m[i] = 0.0
            es.tau_m[i] = 0.0
    es.t_m += dt
    return es


def sdm_decoherence(
    es: ElectronicStructure, dt: float, act_st: int, rates: np.ndarray, tol: float = 0.0
) -> ElectronicStructure:
    if act_st < 0 or act_st >= es.num_states:
        raise ValueError("Error in sdm_decoherence: active state index out of range")
    es = update_populations(es)
    p_aa_old = float(np.real(es.A[act_st, act_st]))
    if p_aa_old > 1.0 + tol:
        sclf = 1.0 / math.sqrt(p_aa_old)
        es.Ccurr *= sclf
        es = update_populations(es)
        p_aa_old = float(np.real(es.A[act_st, act_st]))
    if p_aa_old <= 0.0:
        return es

    inact_st_pop = 0.0
    for i in range(es.num_states):
        if i != act_st:
            itau = float(np.real(rates[i, act_st]))
            sclf = math.exp(-dt * itau)
            es.Ccurr[i] *= sclf
            inact_st_pop += abs(es.Ccurr[i]) ** 2
    if inact_st_pop > 1.0:
        raise ValueError("Error in sdm_decoherence: inactive-state population > 1.0")
    p_aa_new = 1.0 - inact_st_pop
    if p_aa_new < 0.0:
        raise ValueError(
            "Error in sdm_decoherence: new active-state population is negative"
        )
    sclf = math.sqrt(p_aa_new / p_aa_old)
    es.Ccurr[act_st] *= sclf
    es = update_populations(es)
    new_norm = float(np.real(np.vdot(es.Ccurr, es.Ccurr)))
    if abs(new_norm - 1.0) > 0.1:
        raise ValueError("Error in sdm_decoherence: norm deviates too much from 1")
    return es


def gfsh_decohere(
    es: ElectronicStructure, dt: float, temp: float
) -> ElectronicStructure:
    es = update_populations(es)
    curr_state = es.curr_state
    pop_diag = np.real(np.diag(es.A))
    total_pop = np.sum(pop_diag)
    if total_pop <= 0:
        return es
    normalized_pop = pop_diag / total_pop
    num_decohere = int(np.sum(normalized_pop > np.random.rand(*normalized_pop.shape)))
    if num_decohere <= 0:
        es.t_m += dt
        return es
    decohere_states = np.random.permutation(es.num_states)[
        : min(num_decohere, es.num_states)
    ]
    for i in decohere_states:
        i = int(i)
        zeta = np.random.rand()
        P = float(np.real(es.A[i, i])) / total_pop
        dE = float(np.real(es.Hcurr[i, i] - es.Hcurr[curr_state, curr_state]))
        if dE > 0:
            P *= math.exp(-(dE / (kb * temp)))
        if zeta < P:
            es = decohere(es, i)
            curr_state = i
        else:
            es = project_out(es, i)
        es.t_m[i] = 0.0
        es.tau_m[i] = 0.0
    es.t_m += dt
    return es


def decoherence_rates(
    x: np.ndarray, dt: float, icond: int, i1: int, j1: int
) -> Tuple[float, np.ndarray]:
    x = np.asarray(x, dtype=float)
    length = len(x)
    sz = length // 2
    C = np.zeros(sz, dtype=float)
    IC = np.zeros(sz, dtype=float)
    IIC = np.zeros(sz, dtype=float)
    T = np.zeros(sz, dtype=float)
    selIIC = np.zeros(sz, dtype=float)

    for t in range(sz):
        C[t] = np.sum(x[:sz] * x[t : sz + t]) / sz

    sum0 = 0.0
    for t in range(sz):
        IC[t] = sum0
        sum0 += C[t] * (dt / hbar)

    sum0 = 0.0
    for t in range(sz):
        IIC[t] = sum0
        sum0 += IC[t] * (dt / hbar)

    D = np.exp(-IIC)
    nrm = C[0] if sz > 0 else 0.0
    if abs(nrm) > np.finfo(float).eps:
        C /= nrm

    dE = 0.0025
    Npoints = 2000
    J = np.zeros(Npoints, dtype=float)
    tt = np.arange(1, sz, dtype=float) * dt
    Ct = C[1:]

    spectral_path = OUTDIR / f"icond{icond}pair{i1}_{j1}Spectral_density.txt"
    with spectral_path.open("w", encoding="utf-8") as out1:
        for w in range(Npoints):
            ww = w * dE
            J[w] = (1.0 + 2.0 * np.sum(np.cos(ww * tt) * Ct)) * dt
            J[w] = J[w] * J[w] / (2.0 * math.pi)
            out1.write(
                f"w(eV)= {w * dE:.12g} w(cm^-1)= {w * dE * 8065.54468111324:.12g} "
                f"J= {J[w]:.12g} sqrt(J)= {math.sqrt(max(J[w], 0.0)):.12g}\n"
            )

    first = True
    cnt = 0
    for t in range(sz):
        if first:
            if IIC[t] < 2.3:
                T[cnt] = t * t * dt * dt
                selIIC[cnt] = IIC[t]
                cnt += 1
            else:
                first = False

    T = T[:cnt]
    selIIC = selIIC[:cnt]
    a = 0.0
    b = 0.0
    if T.size > 0:
        _, _, a, b = regression(T, selIIC, a, b)
    if b < 0.0:
        b = 0.0

    dephasing_path = OUTDIR / f"icond{icond}pair{i1}_{j1}Dephasing_function.txt"
    with dephasing_path.open("w", encoding="utf-8") as out:
        out.write(
            "Time    D(t)       fitted D(t)     Normalized_autocorrelation_function  "
            "Unnormalized_autocorrelation_function   Second cumulant\n"
        )
        for t in range(sz):
            fitted = math.exp(-a) * math.exp(-b * t * t * dt * dt)
            out.write(
                f"{t * dt:.12g}  {D[t]:.12g}  {fitted:.12g}  {C[t]:.12g} {nrm * C[t]:.12g}  {IIC[t]:.12g}\n"
            )

    return math.sqrt(b), x


def run_decoherence_rates(
    me_es: List[ElectronicStructure], icond: int
) -> List[ElectronicStructure]:
    sz = len(me_es)
    N = me_es[0].num_states
    rij = np.zeros((N, N), dtype=float)
    out_path = OUTDIR / f"decoherence_rates_icond{icond}.txt"
    with out_path.open("w", encoding="utf-8") as out:
        for i1 in range(N):
            for j1 in range(N):
                if i1 == j1:
                    rij[i1, j1] = 0.0
                else:
                    Eij = np.zeros(sz, dtype=float)
                    ave_dEij = 0.0
                    for t in range(sz):
                        dEij = float(
                            np.real(me_es[t].Hcurr[i1, i1] - me_es[t].Hcurr[j1, j1])
                        )
                        Eij[t] = dEij
                        ave_dEij += dEij
                    Eij -= ave_dEij / sz
                    rij[i1, j1], _ = decoherence_rates(Eij, nucl_dt, icond, i1, j1)
                out.write(f"{rij[i1, j1]:.12g} ")
            out.write("\n")
    return me_es


def run_namd(
    me_es: List[ElectronicStructure], me_states: List[MEState], icond: int
) -> Tuple[List[ElectronicStructure], List[MEState]]:
    sz = len(me_es)
    nst = me_es[0].num_states
    init_state = me_es[0].curr_state
    sh_pops = np.zeros((sz, nst), dtype=float)
    se_pops = np.zeros((sz, nst), dtype=float)
    rates = np.zeros((nst, nst), dtype=np.complex128)

    if decoherence_method > 0:
        filename = OUTDIR / f"decoherence_rates_icond{icond}.txt"
        r_ij = np.loadtxt(filename)
        rates = r_ij.astype(np.complex128)

    for _ in range(num_sh_traj):
        me_es[0] = set_state(me_es[0], init_state)
        me_es[0].t_m[0] = 0.0
        for i in range(sz):
            if i > 0:
                me_es[i] = insertion(me_es[i], me_es[i - 1])
            me_es[i] = init_hop_prob1(me_es[i])
            me_es, rates = propagate_electronic(me_es, i, rates)
            me_es[i] = update_populations(me_es[i])

            if decoherence_method == 0:
                me_es[i].g, me_es[i].curr_state = hop(
                    me_es[i].g, me_es[i].curr_state, nst
                )
            elif decoherence_method == 1:
                me_es[i] = dish_decoherence(me_es[i], nucl_dt, 1, Temp, rates)
            elif decoherence_method == 2:
                me_es[i].g, me_es[i].curr_state = hop(
                    me_es[i].g, me_es[i].curr_state, nst
                )
                curr_state = me_es[i].curr_state
                me_es[i] = sdm_decoherence(me_es[i], nucl_dt, curr_state, rates)
                me_es[i] = update_populations(me_es[i])
            elif decoherence_method == 3:
                me_es[i] = gfsh_decohere(me_es[i], nucl_dt, Temp)

            curr_state = me_es[i].curr_state
            sh_pops[i, curr_state] += 1
            se_pops[i, :] += np.real(np.diag(me_es[i].A))

    se_pops2 = se_pops / num_sh_traj
    sh_pops2 = sh_pops / num_sh_traj

    outfile1 = OUTDIR / f"me_pop{icond}"
    outfile2 = OUTDIR / f"out{icond}"
    with outfile1.open("w", encoding="utf-8") as out1, outfile2.open(
        "w", encoding="utf-8"
    ) as out2:
        for i in range(sz):
            out1.write(f"time {i} ")
            tot = np.sum(se_pops2[i, :])
            for j in range(nst):
                out1.write(f"P({j})= {se_pops2[i, j]:.10f} ")
            out1.write(f"Total= {tot:.10f}\n")

            out2.write(f"time {i} ")
            for j in range(nst):
                out2.write(f"P({j})= {sh_pops2[i, j]:.10f} ")
            out2.write("\n")

    return me_es, me_states


def load_hamiltonian_batch() -> List[np.ndarray]:
    H_batch: List[np.ndarray] = []
    for j in range(Ham_size):
        ham_re_path = RT / "res" / f"0_Ham_{j}_re"
        ham_im_path = RT / "res" / f"0_Ham_{j}_im"
        Ham_re = np.loadtxt(ham_re_path)
        Ham_im = np.loadtxt(ham_im_path)
        Ham = extract_2d(Ham_re, active_space, -1) + 1j * extract_2d(
            Ham_im, active_space, -1
        )
        H_batch.append(Ham * Ry_to_eV)
    return H_batch


def main() -> None:
    start_time = time.time()

    me_states = [make_me_state(s) for s in states_def]
    me_numstates = len(me_states)
    numstates = len(active_space)
    num_elec = len(me_states[0].actual_state)

    iconds = np.zeros((iconds_i, 2), dtype=int)
    for ii in range(iconds_i):
        iconds[ii, :] = [ii, me_numstates - 1]

    H_batch = load_hamiltonian_batch()

    for icond in range(iconds.shape[0]):
        oe_es = [ElectronicStructure(2 * numstates) for _ in range(namdtime)]
        me_es = [ElectronicStructure(me_numstates) for _ in range(namdtime)]

        istart = int(iconds[icond, 0])
        for j in range(istart, istart + namdtime):
            t = j - istart
            j1 = j % Ham_size
            Hij = H_batch[j1]
            Htmp = np.zeros((2 * numstates, 2 * numstates), dtype=np.complex128)
            for k1 in range(numstates):
                i1 = 2 * k1
                i2 = 2 * k1 + 1
                for k2 in range(numstates):
                    jx1 = 2 * k2
                    jx2 = 2 * k2 + 1
                    Htmp[i1, jx1] = Hij[k1, k2]
                    Htmp[i2, jx2] = Hij[k1, k2]

            if alp_bet != 0:
                for k1 in range(numstates):
                    i1 = 2 * k1
                    i2 = 2 * k1 + 1
                    for k2 in range(numstates):
                        jx1 = 2 * k2
                        jx2 = 2 * k2 + 1
                        Htmp[i2, jx1] = Hij[k1, k2]
                        Htmp[i1, jx2] = Hij[k1, k2]
            oe_es[t].Hcurr = Htmp

        me_es[0] = set_state(me_es[0], int(iconds[icond, 1]))

        for j in range(istart, istart + namdtime):
            t = j - istart
            Hme = np.zeros((me_numstates, me_numstates), dtype=np.complex128)
            for I in range(me_numstates):
                for J in range(me_numstates):
                    orb_i = 0
                    orb_j = 0
                    delt, _, _, orb_i, orb_j = delta(
                        me_states[I].actual_state,
                        me_states[J].actual_state,
                        orb_i,
                        orb_j,
                    )
                    if delt:
                        oi = ext2int(orb_i)
                        oj = ext2int(orb_j)
                        Hme[I, J] += oe_es[t].Hcurr[oi, oj]

                diagE = 0.0 + 0.0j
                for el in range(num_elec):
                    oi = ext2int(int(me_states[I].actual_state[el]))
                    diagE += oe_es[t].Hcurr[oi, oi]
                shift = me_states[I].shift
                Hme[I, I] = diagE + shift
            me_es[t].Hcurr = Hme

        outfile = OUTDIR / f"me_energies{icond}"
        with outfile.open("w", encoding="utf-8") as out:
            for j in range(istart, istart + namdtime):
                t = j - istart
                e0 = me_es[t].Hcurr[0, 0]
                out.write(f"t= {j:d}  E[0]= {np.real(e0):.12g}  ")
                for I in range(me_es[t].num_states):
                    out.write(
                        f"E[{I}]-E[0]= {np.real(me_es[t].Hcurr[I, I] - e0):.12g}  "
                    )
                out.write("\n")

        if decoherence_method > 0:
            me_es = run_decoherence_rates(me_es, icond)

        me_es, me_states = run_namd(me_es, me_states, icond)

        del oe_es
        del me_es

    elapsed_time = time.time() - start_time
    with (SCRIPT_DIR / "log.txt").open("w", encoding="utf-8") as log_file:
        print(f"Running time: {elapsed_time:.1f} seconds")
        log_file.write(f"Running time: {elapsed_time:.1f} seconds\n")

    me_pop0 = OUTDIR / "me_pop0"
    last_line = ""
    if me_pop0.exists():
        with me_pop0.open("r", encoding="utf-8") as f:
            for line in f:
                last_line = line.rstrip("\n")
    print(last_line)


if __name__ == "__main__":
    main()
