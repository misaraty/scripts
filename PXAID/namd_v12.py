from __future__ import annotations

import gc
import math
import random
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


RANDOM_SEED = 42
Ry_to_eV = 13.60569253
Ha_to_eV = 27.211396

random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)


@dataclass
class NAMDParams:
    hbar: float
    kb: float
    nucl_dt: float
    elec_dt: float
    integrator: int
    sh_algo: int
    num_sh_traj: int
    decoherence_method: int
    boltz_flag: int
    Temp: float


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
        n = int(self.num_states)
        self.num_states = n
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
    active_space: List[int]
    actual_state: List[int]
    shift: float
    nac_scl_indx: List[int] = field(default_factory=list)
    nac_scl: List[float] = field(default_factory=list)


def make_me_state(cs: Sequence, active_space: Sequence[int]) -> MEState:
    shift = float(cs[2]) if len(cs) >= 3 else 0.0
    return MEState(str(cs[0]), list(active_space), list(cs[1]), shift)


def load_matrix(path: str | Path) -> np.ndarray:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Cannot find file: {path}")
    data = np.loadtxt(path)
    return np.asarray(data, dtype=float)


def extract_2d(
    input_matrix: np.ndarray, templ: Sequence[int], shift: int
) -> np.ndarray:
    idx = np.asarray(templ, dtype=int) + int(shift)
    return input_matrix[np.ix_(idx, idx)]


def build_active_map(active_space: Sequence[int]) -> Dict[int, int]:
    return {int(el): k for k, el in enumerate(active_space)}


def ext2int(external: int, active_map: Dict[int, int]) -> int:
    idx = active_map[abs(int(external))]
    f = 1 if external < 0 else 0
    return 2 * abs(idx) + f


def delta_states(
    A: Sequence[int], B: Sequence[int], a: int = 0, b: int = 0
) -> Tuple[bool, Sequence[int], Sequence[int], int, int]:
    C = set(A) | set(B)
    nexc = 0
    for c in C:
        n_in_a = sum(1 for x in A if x == c)
        n_in_b = sum(1 for x in B if x == c)
        d = n_in_a - n_in_b
        if d > 0:
            nexc += d
        if d == 1:
            a = c
        if d == -1:
            b = c
    return nexc == 1, A, B, a, b


def regression(
    X: np.ndarray, Y: np.ndarray, a: float = 0.0, b: float = 0.0
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    sx = float(np.sum(X))
    if abs(sx) < np.finfo(float).eps:
        b = 0.0
    else:
        b = float(np.sum(Y) / sx)
    return X, Y, a, b


def set_state(es: ElectronicStructure, indx: int) -> ElectronicStructure:
    es.Ccurr.fill(0.0 + 0.0j)
    es.Ccurr[indx] = 1.0 + 0.0j
    es.curr_state = int(indx)
    return es


def init_hop_prob1(es: ElectronicStructure) -> ElectronicStructure:
    es.g.fill(0.0)
    np.fill_diagonal(es.g, 1.0)
    return es


def efield(t: float, E: np.ndarray, Eex: float) -> Tuple[np.ndarray, float]:
    E.fill(0.0)
    return E, 0.0


def build_transition_map(
    me_states: Sequence[MEState], active_map: Dict[int, int]
) -> List[Tuple[int, int, int, int]]:
    me_numstates = len(me_states)
    transition_map: List[Tuple[int, int, int, int]] = []

    for Istate in range(me_numstates):
        for Jstate in range(me_numstates):
            orb_i = 0
            orb_j = 0
            delt, _, _, orb_i, orb_j = delta_states(
                me_states[Istate].actual_state,
                me_states[Jstate].actual_state,
                orb_i,
                orb_j,
            )
            if delt:
                ii = ext2int(orb_i, active_map)
                jj = ext2int(orb_j, active_map)
                transition_map.append((Istate, Jstate, ii, jj))

    return transition_map


def build_diag_orbital_map(
    me_states: Sequence[MEState], active_map: Dict[int, int]
) -> List[List[int]]:
    return [[ext2int(el, active_map) for el in st.actual_state] for st in me_states]


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
    es: ElectronicStructure, Hij: complex, dt: float, i: int, j: int, p: NAMDParams
) -> ElectronicStructure:
    phi1 = 0.5 * dt * Hij.imag / p.hbar
    phi2 = -dt * Hij.real / p.hbar
    rot1(es, phi1, i, j)
    rot2(es, phi2, i, j)
    rot1(es, phi1, i, j)
    return es


def phase(
    es: ElectronicStructure, Hii: complex, dt: float, i: int, p: NAMDParams
) -> ElectronicStructure:
    phi = -dt * Hii.real / p.hbar
    es.Ccurr[i] = np.exp(1j * phi) * es.Ccurr[i]
    return es


def propagate_coefficients(
    es: ElectronicStructure, dt: float, Ef: np.ndarray, p: NAMDParams
) -> ElectronicStructure:
    for i in range(es.num_states):
        for j in range(i + 1, es.num_states):
            rot(es, es.Hcurr[i, j], 0.5 * dt, i, j, p)

    for i in range(es.num_states):
        phase(es, es.Hcurr[i, i], dt, i, p)

    for i in range(es.num_states - 1, -1, -1):
        for j in range(es.num_states - 1, i, -1):
            rot(es, es.Hcurr[i, j], 0.5 * dt, i, j, p)

    return es


def update_populations(es: ElectronicStructure) -> ElectronicStructure:
    es.A[:, :] = np.outer(np.conjugate(es.Ccurr), es.Ccurr)
    return es


def update_hop_prob_fssh(
    es: ElectronicStructure,
    dt: float,
    boltz_flag_arg: int,
    Temp_arg: float,
    Ef: np.ndarray,
    Eex: float,
    rates: np.ndarray,
    p: NAMDParams,
) -> ElectronicStructure:
    update_populations(es)
    Heff = es.Hcurr

    for i in range(es.num_states):
        a_ii = float(es.A[i, i].real)
        if a_ii < 1e-12:
            a_ii = 1e-12

        g_row_sum = 0.0
        for j in range(es.num_states):
            if j != i:
                g_ij = (2.0 * dt / (a_ii * p.hbar)) * (es.A[i, j] * Heff[i, j]).imag
                if g_ij < 0.0:
                    g_ij = 0.0

                E_i = Heff[i, i].real
                E_j = Heff[j, j].real
                dE = E_j - E_i
                bf = 1.0

                if boltz_flag_arg != 0 and dE > Eex:
                    bf = math.exp(-((dE - Eex) / (p.kb * Temp_arg)))

                g_ij *= bf
                g_row_sum += g_ij
                es.g[i, j] = g_ij
        es.g[i, i] = 1.0 - g_row_sum

    return es


def propagate_electronic(
    es: ElectronicStructure, istep: int, rates: np.ndarray, p: NAMDParams
) -> ElectronicStructure:
    nel = int(round(p.nucl_dt / p.elec_dt))
    Eex = 0.0
    Ef = np.zeros(3, dtype=float)

    if p.integrator == 0:
        for j in range(nel):
            tim = istep * p.nucl_dt + j * p.elec_dt
            Ef, Eex = efield(tim, Ef, Eex)
            propagate_coefficients(es, p.elec_dt, Ef, p)
            if p.sh_algo == 0:
                update_hop_prob_fssh(
                    es, p.elec_dt, p.boltz_flag, p.Temp, Ef, Eex, rates, p
                )

    return es


def hop(sh_prob: np.ndarray, hopstate: int) -> int:
    input_state = int(hopstate)
    ksi = random.random()
    probs = sh_prob[input_state, :]
    nrm = float(np.sum(probs))
    if nrm <= 0:
        return input_state

    accum = 0.0
    for j, prob in enumerate(probs):
        accum += float(prob) / nrm
        if accum > 0 and ksi <= accum:
            return j

    return input_state


def nextpow2_int(n: int) -> int:
    return 1 << int(math.ceil(math.log2(n)))


def autocorr_fft_cross_strict(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    length = len(x)
    sz = int(math.floor(length / 2))

    if sz <= 0:
        return np.array([], dtype=float)
    if length < 2 * sz:
        raise ValueError(
            "Input length is inconsistent with strict autocorrelation definition."
        )

    a = x[:sz]
    b = x[: 2 * sz - 1]

    nconv = len(a) + len(b) - 1
    nfft = nextpow2_int(nconv)

    ar = a[::-1].copy()
    avec = np.concatenate([ar, np.zeros(nfft - len(ar), dtype=float)])
    bvec = np.concatenate([b.copy(), np.zeros(nfft - len(b), dtype=float)])

    conv_full = np.fft.ifft(np.fft.fft(avec) * np.fft.fft(bvec)).real

    C = np.zeros(sz, dtype=float)
    start_idx = sz - 1
    for t in range(sz):
        C[t] = conv_full[start_idx + t] / sz

    return C


def decoherence_rates(
    x: np.ndarray,
    dt: float,
    icond: int,
    i1: int,
    j1: int,
    outdir: str | Path,
    p: NAMDParams,
) -> Tuple[float, np.ndarray]:
    x = np.asarray(x, dtype=float)
    sz = int(math.floor(len(x) / 2))
    outdir = Path(outdir)

    C = autocorr_fft_cross_strict(x)
    IC = np.zeros(sz, dtype=float)
    IIC = np.zeros(sz, dtype=float)
    T = np.zeros(sz, dtype=float)
    selIIC = np.zeros(sz, dtype=float)

    sum0 = 0.0
    for t in range(sz):
        IC[t] = sum0
        sum0 += C[t] * (dt / p.hbar)

    sum0 = 0.0
    for t in range(sz):
        IIC[t] = sum0
        sum0 += IC[t] * (dt / p.hbar)

    D = np.exp(-IIC)

    nrm = C[0]
    if abs(nrm) > np.finfo(float).eps:
        C = C / nrm

    dE = 0.0025
    Npoints = 2000
    J = np.zeros(Npoints, dtype=float)
    tt = np.arange(1, sz, dtype=float) * dt
    Ct = C[1:]

    spectral_path = outdir / f"icond{icond}pair{i1}_{j1}Spectral_density.txt"
    with spectral_path.open("w", encoding="utf-8") as io:
        for w in range(Npoints):
            ww = w * dE
            J[w] = (1.0 + 2.0 * np.sum(np.cos(ww * tt) * Ct)) * dt
            J[w] = J[w] ** 2 / (2.0 * math.pi)
            io.write(
                f"w(eV)= {w * dE:.12g} w(cm^-1)= {w * dE * 8065.54468111324:.12g} "
                f"J= {J[w]:.12g} sqrt(J)= {math.sqrt(J[w]):.12g}\n"
            )

    first_region = True
    cnt = 0
    for t in range(sz):
        if first_region:
            if IIC[t] < 2.3:
                T[cnt] = t * t * dt * dt
                selIIC[cnt] = IIC[t]
                cnt += 1
            else:
                first_region = False

    Tfit = T[:cnt]
    selIICfit = selIIC[:cnt]

    a = 0.0
    b = 0.0
    if cnt > 0:
        _, _, a, b = regression(Tfit, selIICfit, a, b)
    if b < 0.0:
        b = 0.0

    dephase_path = outdir / f"icond{icond}pair{i1}_{j1}Dephasing_function.txt"
    with dephase_path.open("w", encoding="utf-8") as io:
        io.write(
            "Time    D(t)       fitted D(t)     Normalized_autocorrelation_function  "
            "Unnormalized_autocorrelation_function   Second cumulant\n"
        )
        for t in range(sz):
            io.write(
                f"{t * dt:.12g}  {D[t]:.12g}  {math.exp(-a) * math.exp(-b * t * t * dt * dt):.12g}  "
                f"{C[t]:.12g} {nrm * C[t]:.12g}  {IIC[t]:.12g}\n"
            )

    return math.sqrt(b), x


def run_decoherence_rates(
    Hme_batch: Sequence[np.ndarray], icond: int, outdir: str | Path, p: NAMDParams
) -> np.ndarray:
    sz = len(Hme_batch)
    N = Hme_batch[0].shape[0]
    rij = np.zeros((N, N), dtype=float)
    outdir = Path(outdir)
    deco_path = outdir / f"decoherence_rates_icond{icond}.txt"

    with deco_path.open("w", encoding="utf-8") as io:
        for i1 in range(N):
            for j1 in range(N):
                if i1 == j1:
                    rij[i1, j1] = 0.0
                else:
                    Eij = np.zeros(sz, dtype=float)
                    ave_dEij = 0.0

                    for t in range(sz):
                        dEij = (Hme_batch[t][i1, i1] - Hme_batch[t][j1, j1]).real
                        Eij[t] = dEij
                        ave_dEij += dEij

                    Eij -= ave_dEij / sz
                    rij[i1, j1], _ = decoherence_rates(
                        Eij, p.nucl_dt, icond, i1, j1, outdir, p
                    )
                io.write(f"{rij[i1, j1].real:.12g} ")
            io.write("\n")

    return rij


def update_decoherence_times(
    es: ElectronicStructure, rates: np.ndarray
) -> ElectronicStructure:
    update_populations(es)
    pops = np.real(np.diag(es.A))
    es.tau_m[:] = np.real(rates) @ pops
    return es


def decohere(es: ElectronicStructure, i: int) -> ElectronicStructure:
    es.Ccurr.fill(0.0 + 0.0j)
    es.Ccurr[i] = 1.0 + 0.0j
    es.curr_state = int(i)
    return es


def project_out(es: ElectronicStructure, i: int) -> ElectronicStructure:
    update_populations(es)
    es.Ccurr[i] = 0.0 + 0.0j
    nrm2 = float(np.sum(np.real(np.diag(es.A))) - es.A[i, i].real)
    if nrm2 <= 0:
        es.Ccurr.fill(0.0 + 0.0j)
        es.Ccurr[es.curr_state] = 1.0 + 0.0j
    else:
        es.Ccurr /= math.sqrt(nrm2)
    update_populations(es)
    return es


def dish_decoherence(
    es: ElectronicStructure,
    dt: float,
    boltz_flag_arg: int,
    Temp_arg: float,
    rates: np.ndarray,
    p: NAMDParams,
) -> ElectronicStructure:
    update_decoherence_times(es, rates)

    for i in range(es.num_states):
        if es.tau_m[i] <= 0:
            continue

        rnd_i = 1.0 / es.tau_m[i]
        if es.t_m[i] >= rnd_i:
            zeta = random.random()
            P = es.A[i, i].real
            dE = (es.Hcurr[i, i] - es.Hcurr[es.curr_state, es.curr_state]).real
            if dE > 0:
                P *= math.exp(-(dE / (p.kb * Temp_arg)))

            if zeta < P:
                decohere(es, i)
                break
            else:
                project_out(es, i)
                es.t_m[i] = 0.0
                es.tau_m[i] = 0.0

    es.t_m += dt
    return es


def sdm_decoherence(
    es: ElectronicStructure,
    dt: float,
    act_st: int,
    rates: np.ndarray,
    tol: float = 0.0,
) -> ElectronicStructure:
    if act_st < 0 or act_st >= es.num_states:
        raise IndexError(
            "Error in ElectronicStructure_sdm_decoherence: active state index out of range"
        )

    update_populations(es)
    p_aa_old = es.A[act_st, act_st].real

    if p_aa_old > 1.0 + tol:
        sclf = 1.0 / math.sqrt(p_aa_old)
        es.Ccurr *= sclf
        update_populations(es)
        p_aa_old = es.A[act_st, act_st].real

    if p_aa_old <= 0.0:
        return es

    inact_st_pop = 0.0
    for i in range(es.num_states):
        if i != act_st:
            itau = rates[i, act_st].real
            sclf = math.exp(-dt * itau)
            es.Ccurr[i] *= sclf
            inact_st_pop += abs(es.Ccurr[i]) ** 2

    if inact_st_pop > 1.0:
        raise ValueError(
            "Error in ElectronicStructure_sdm_decoherence: inactive-state population > 1.0"
        )

    p_aa_new = 1.0 - inact_st_pop
    if p_aa_new < 0.0:
        raise ValueError(
            "Error in ElectronicStructure_sdm_decoherence: new active-state population is negative"
        )

    sclf = math.sqrt(p_aa_new / p_aa_old)
    es.Ccurr[act_st] *= sclf
    update_populations(es)

    new_norm = float(np.real(np.vdot(es.Ccurr, es.Ccurr)))
    if abs(new_norm - 1.0) > 0.1:
        raise ValueError(
            "Error in ElectronicStructure_sdm_decoherence: norm deviates too much from 1"
        )

    return es


def gfsh_decohere(
    es: ElectronicStructure, dt: float, Temp_arg: float, p: NAMDParams
) -> ElectronicStructure:
    update_populations(es)
    curr_state = es.curr_state
    pop_diag = np.real(np.diag(es.A))
    total_pop = float(np.sum(pop_diag))

    if total_pop <= 0:
        return es

    normalized_pop = pop_diag / total_pop
    num_decohere = int(np.sum(normalized_pop > np.random.random(len(normalized_pop))))

    if num_decohere <= 0:
        es.t_m += dt
        return es

    decohere_states = (
        np.random.permutation(es.num_states)[: min(num_decohere, es.num_states)]
    ).tolist()

    for i in decohere_states:
        zeta = random.random()
        P = es.A[i, i].real / total_pop
        dE = (es.Hcurr[i, i] - es.Hcurr[curr_state, curr_state]).real

        if dE > 0:
            P *= math.exp(-(dE / (p.kb * Temp_arg)))

        if zeta < P:
            decohere(es, i)
            curr_state = i
        else:
            project_out(es, i)

        es.t_m[i] = 0.0
        es.tau_m[i] = 0.0

    es.t_m += dt
    return es


def build_spin_orbital_H(Htmp: np.ndarray, Hij: np.ndarray, alp_bet: int) -> np.ndarray:
    numstates = Hij.shape[0]
    Htmp.fill(0.0 + 0.0j)

    for k1 in range(numstates):
        i1 = 2 * k1
        i2 = 2 * k1 + 1
        for k2 in range(numstates):
            jx1 = 2 * k2
            jx2 = 2 * k2 + 1
            val = Hij[k1, k2]
            Htmp[i1, jx1] = val
            Htmp[i2, jx2] = val

    if alp_bet != 0:
        for k1 in range(numstates):
            i1 = 2 * k1
            i2 = 2 * k1 + 1
            for k2 in range(numstates):
                jx1 = 2 * k2
                jx2 = 2 * k2 + 1
                val = Hij[k1, k2]
                Htmp[i2, jx1] = val
                Htmp[i1, jx2] = val

    return Htmp


def build_Hme_batch(
    H_batch: Sequence[np.ndarray],
    istart: int,
    namdtime: int,
    Ham_size: int,
    me_numstates: int,
    me_states: Sequence[MEState],
    transition_map: Sequence[Tuple[int, int, int, int]],
    diag_orbital_map: Sequence[Sequence[int]],
    alp_bet: int,
) -> List[np.ndarray]:
    numstates = H_batch[0].shape[0]
    Htmp = np.zeros((2 * numstates, 2 * numstates), dtype=np.complex128)
    Hme = np.zeros((me_numstates, me_numstates), dtype=np.complex128)
    Hme_batch: List[np.ndarray] = [None] * namdtime  # type: ignore[list-item]

    for j in range(istart, istart + namdtime):
        t = j - istart
        j1 = j % Ham_size
        Hij = H_batch[j1]

        build_spin_orbital_H(Htmp, Hij, alp_bet)
        Hme.fill(0.0 + 0.0j)

        for I, J, ii, jj in transition_map:
            Hme[I, J] += Htmp[ii, jj]

        for Istate in range(me_numstates):
            diagE = 0.0 + 0.0j
            for orb in diag_orbital_map[Istate]:
                diagE += Htmp[orb, orb]
            Hme[Istate, Istate] = diagE + me_states[Istate].shift

        Hme_batch[t] = Hme.copy()

    return Hme_batch


def write_me_energies(
    Hme_batch: Sequence[np.ndarray], istart: int, icond: int, outdir: str | Path
) -> None:
    outdir = Path(outdir)
    outfile = outdir / f"me_energies{icond}"
    sz = len(Hme_batch)
    nst = Hme_batch[0].shape[0]

    with outfile.open("w", encoding="utf-8") as io:
        for local_t in range(sz):
            j = istart + local_t
            H = Hme_batch[local_t]
            io.write(f"t= {j:d}  E[0]= {H[0, 0].real:.12g}  ")
            e0 = H[0, 0]
            for Istate in range(nst):
                io.write(f"E[{Istate:d}]-E[0]= {(H[Istate, Istate] - e0).real:.12g}  ")
            io.write("\n")


def run_namd(
    Hme_batch: Sequence[np.ndarray],
    me_states: Sequence[MEState],
    icond: int,
    init_state: int,
    outdir: str | Path,
    p: NAMDParams,
) -> None:
    sz = len(Hme_batch)
    nst = Hme_batch[0].shape[0]
    sh_pops = np.zeros((sz, nst), dtype=float)
    se_pops = np.zeros((sz, nst), dtype=float)
    rates = np.zeros((nst, nst), dtype=np.complex128)
    outdir = Path(outdir)

    if p.decoherence_method > 0:
        filename = outdir / f"decoherence_rates_icond{icond}.txt"
        r_ij = load_matrix(filename)
        rates[:, :] = r_ij.astype(np.complex128)

    es = ElectronicStructure(nst)

    for _ in range(p.num_sh_traj):
        set_state(es, init_state)
        es.t_m.fill(0.0)
        es.tau_m.fill(0.0)

        for i in range(sz):
            es.Hcurr[:, :] = Hme_batch[i]
            init_hop_prob1(es)
            propagate_electronic(es, i, rates, p)
            update_populations(es)

            if p.decoherence_method == 0:
                es.curr_state = hop(es.g, es.curr_state)
            elif p.decoherence_method == 1:
                dish_decoherence(es, p.nucl_dt, 1, p.Temp, rates, p)
            elif p.decoherence_method == 2:
                es.curr_state = hop(es.g, es.curr_state)
                sdm_decoherence(es, p.nucl_dt, es.curr_state, rates)
                update_populations(es)
            elif p.decoherence_method == 3:
                gfsh_decohere(es, p.nucl_dt, p.Temp, p)

            curr_state = es.curr_state
            sh_pops[i, curr_state] += 1.0
            for j in range(nst):
                se_pops[i, j] += es.A[j, j].real

    se_pops2 = se_pops / p.num_sh_traj
    sh_pops2 = sh_pops / p.num_sh_traj
    outfile1 = outdir / f"me_pop{icond}"
    outfile2 = outdir / f"out{icond}"

    with outfile1.open("w", encoding="utf-8") as out1, outfile2.open(
        "w", encoding="utf-8"
    ) as out2:
        for i in range(sz):
            out1.write(f"time {i:d} ")
            tot = float(np.sum(se_pops2[i, :]))
            for j in range(nst):
                out1.write(f"P({j:d})= {se_pops2[i, j]:.10f} ")
            out1.write(f"Total= {tot:.10f}\n")

            out2.write(f"time {i:d} ")
            for j in range(nst):
                out2.write(f"P({j:d})= {sh_pops2[i, j]:.10f} ")
            out2.write("\n")


def main() -> None:
    t0 = time.time()

    filepath = Path(__file__).resolve().parent
    rt = filepath.parent
    outdir = filepath / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    namdtime = 10000
    num_sh_traj = 1000
    iconds_i = 1000
    Ham_size = 4001

    active_space = [1, 2]
    states = [("GS", [1, -1], 0.00), ("S1", [1, -2], 0.00)]

    decoherence_method = 2
    nucl_dt = 1.0
    elec_dt = 1.0
    integrator = 0
    sh_algo = 0
    boltz_flag = 1
    Temp = 300.0
    alp_bet = 0

    p = NAMDParams(
        0.658218,
        8.617e-5,
        nucl_dt,
        elec_dt,
        integrator,
        sh_algo,
        num_sh_traj,
        decoherence_method,
        boltz_flag,
        Temp,
    )

    me_numstates = len(states)
    me_states = [make_me_state(s, active_space) for s in states]

    active_map = build_active_map(active_space)
    transition_map = build_transition_map(me_states, active_map)
    diag_orbital_map = build_diag_orbital_map(me_states, active_map)

    iconds = np.zeros((iconds_i, 2), dtype=int)
    for ii in range(iconds_i):
        iconds[ii, :] = [ii, me_numstates - 1]

    H_batch: List[np.ndarray] = [None] * Ham_size  # type: ignore[list-item]
    for j in range(Ham_size):
        Ham_re = load_matrix(rt / "res" / f"0_Ham_{j}_re")
        Ham_im = load_matrix(rt / "res" / f"0_Ham_{j}_im")
        Ham = extract_2d(Ham_re, active_space, -1) + 1j * extract_2d(
            Ham_im, active_space, -1
        )
        H_batch[j] = Ham * Ry_to_eV

    for icond in range(iconds.shape[0]):
        istart = int(iconds[icond, 0])
        init_state = int(iconds[icond, 1])

        Hme_batch = build_Hme_batch(
            H_batch,
            istart,
            namdtime,
            Ham_size,
            me_numstates,
            me_states,
            transition_map,
            diag_orbital_map,
            alp_bet,
        )

        write_me_energies(Hme_batch, istart, icond, outdir)

        if p.decoherence_method > 0:
            run_decoherence_rates(Hme_batch, icond, outdir, p)

        run_namd(Hme_batch, me_states, icond, init_state, outdir, p)

        del Hme_batch
        gc.collect()

    elapsed_time = time.time() - t0
    with (filepath / "log.txt").open("w", encoding="utf-8") as log_file:
        print(f"Running time: {elapsed_time:.1f} seconds")
        log_file.write(f"Running time: {elapsed_time:.1f} seconds\n")

    mepop0 = outdir / "me_pop0"
    if mepop0.is_file():
        last_line = ""
        with mepop0.open("r", encoding="utf-8") as io:
            for line in io:
                last_line = line.rstrip("\n")
        print(last_line)


if __name__ == "__main__":
    main()
