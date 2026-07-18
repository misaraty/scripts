using LinearAlgebra
using DelimitedFiles
using Printf
using Random
using Statistics
using FFTW
using Base.Threads

const RANDOM_SEED = 42
const Ry_to_eV = 13.60569253
const Ha_to_eV = 27.211396

const SAVE_DECOHERENCE_RATES = false
const SAVE_ICOND_FILES = false
const SAVE_ME_ENERGIES = false
const SAVE_ME_POP = false

Random.seed!(RANDOM_SEED)

struct NAMDParams
    hbar::Float64
    kb::Float64
    nucl_dt::Float64
    elec_dt::Float64
    integrator::Int
    sh_algo::Int
    num_sh_traj::Int
    decoherence_method::Int
    boltz_flag::Int
    Temp::Float64
end

mutable struct MEState
    name::String
    active_space::Vector{Int}
    actual_state::Vector{Int}
    shift::Float64
    nac_scl_indx::Vector{Int}
    nac_scl::Vector{Float64}
end

function make_me_state(cs, active_space::Vector{Int})
    shift = length(cs) >= 3 ? Float64(cs[3]) : 0.0
    return MEState(cs[1], copy(active_space), copy(cs[2]), shift, Int[], Float64[])
end

function load_matrix(path::AbstractString)
    isfile(path) || error("Cannot find file: $path")
    return Matrix{Float64}(readdlm(path))
end

@inline function extract_2D(input_matrix::AbstractMatrix, templ::Vector{Int}, shift::Int)
    idx = templ .+ shift .+ 1
    return input_matrix[idx, idx]
end

function build_active_map(active_space::Vector{Int})
    active_map = Dict{Int, Int}()
    for k in eachindex(active_space)
        active_map[active_space[k]] = k - 1
    end
    return active_map
end

@inline function ext2int(external::Int, active_map::Dict{Int, Int})
    idx = active_map[abs(external)]
    f = external < 0 ? 1 : 0
    return 2 * abs(idx) + f
end

function delta_states(A::Vector{Int}, B::Vector{Int}, a::Int = 0, b::Int = 0)
    C = unique(vcat(A, B))
    nexc = 0
    for c in C
        n_in_a = count(==(c), A)
        n_in_b = count(==(c), B)
        d = n_in_a - n_in_b
        if d > 0
            nexc += d
        end
        if d == 1
            a = c
        end
        if d == -1
            b = c
        end
    end
    return nexc == 1, A, B, a, b
end

function regression(
    X::AbstractVector,
    Y::AbstractVector,
    a::Float64 = 0.0,
    b::Float64 = 0.0,
)
    sx = sum(X)
    b = abs(sx) < eps(Float64) ? 0.0 : sum(Y) / sx
    return X, Y, a, b
end

function build_transition_map(me_states::Vector{MEState}, active_map::Dict{Int, Int})
    me_numstates = length(me_states)
    transition_map = Vector{NTuple{4, Int}}()
    for Istate in 0:(me_numstates - 1), Jstate in 0:(me_numstates - 1)
        orb_i = 0
        orb_j = 0
        delt, _, _, orb_i, orb_j = delta_states(
            me_states[Istate + 1].actual_state,
            me_states[Jstate + 1].actual_state,
            orb_i,
            orb_j,
        )
        if delt
            ii = ext2int(orb_i, active_map) + 1
            jj = ext2int(orb_j, active_map) + 1
            push!(transition_map, (Istate + 1, Jstate + 1, ii, jj))
        end
    end
    return transition_map
end

function build_diag_orbital_map(me_states::Vector{MEState}, active_map::Dict{Int, Int})
    return [[ext2int(el, active_map) + 1 for el in st.actual_state] for st in me_states]
end

function pack_diag_orbital_map(diag_orbital_map::Vector{Vector{Int}})
    nst = length(diag_orbital_map)
    maxlen = maximum(length.(diag_orbital_map))
    packed = zeros(Int, nst, maxlen)
    @inbounds for i in 1:nst, k in 1:length(diag_orbital_map[i])
        packed[i, k] = diag_orbital_map[i][k]
    end
    return packed
end

@inline function zero_state!(C::Vector{ComplexF64}, idx0::Int)
    fill!(C, 0.0 + 0.0im)
    C[idx0 + 1] = 1.0 + 0.0im
    return nothing
end

@inline function update_populations_core!(C::Vector{ComplexF64}, A::Matrix{ComplexF64})
    n = length(C)
    @inbounds for i in 1:n
        ci = conj(C[i])
        for j in 1:n
            A[i, j] = ci * C[j]
        end
    end
    return nothing
end

@inline function init_hop_prob_core!(g::Matrix{Float64})
    fill!(g, 0.0)
    @inbounds for i in 1:size(g, 1)
        g[i, i] = 1.0
    end
    return nothing
end

@inline function rot1_core!(C::Vector{ComplexF64}, phi::Float64, i::Int, j::Int)
    c = cos(phi)
    s = sin(phi)
    ci = C[i]
    cj = C[j]
    C[i] = c * ci + s * cj
    C[j] = -s * ci + c * cj
    return nothing
end

@inline function rot2_core!(C::Vector{ComplexF64}, phi::Float64, i::Int, j::Int)
    cs = cos(phi)
    isi = 1im * sin(phi)
    ci = C[i]
    cj = C[j]
    C[i] = cs * ci + isi * cj
    C[j] = isi * ci + cs * cj
    return nothing
end

@inline function rot_core!(
    C::Vector{ComplexF64},
    Hij::ComplexF64,
    dt::Float64,
    i::Int,
    j::Int,
    hbar::Float64,
)
    phi1 = 0.5 * dt * imag(Hij) / hbar
    phi2 = -dt * real(Hij) / hbar
    rot1_core!(C, phi1, i, j)
    rot2_core!(C, phi2, i, j)
    rot1_core!(C, phi1, i, j)
    return nothing
end

function propagate_coefficients_core!(
    C::Vector{ComplexF64},
    H::AbstractMatrix{ComplexF64},
    dt::Float64,
    hbar::Float64,
)
    n = length(C)
    halfdt = 0.5 * dt
    @inbounds for i in 1:n
        for j in (i + 1):n
            rot_core!(C, H[i, j], halfdt, i, j, hbar)
        end
    end
    @inbounds for i in 1:n
        phi = -dt * real(H[i, i]) / hbar
        C[i] = exp(1im * phi) * C[i]
    end
    @inbounds for i in n:-1:1
        for j in n:-1:(i + 1)
            rot_core!(C, H[i, j], halfdt, i, j, hbar)
        end
    end
    return nothing
end

function update_hop_prob_fssh_core!(
    C::Vector{ComplexF64},
    H::AbstractMatrix{ComplexF64},
    A::Matrix{ComplexF64},
    g::Matrix{Float64},
    dt::Float64,
    hbar::Float64,
    kb::Float64,
    boltz_flag::Int,
    temp::Float64,
    eex::Float64,
)
    update_populations_core!(C, A)
    n = length(C)
    @inbounds for i in 1:n
        a_ii = real(A[i, i])
        if a_ii < 1e-12
            a_ii = 1e-12
        end
        g_row_sum = 0.0
        for j in 1:n
            if j != i
                gij = (2.0 * dt / (a_ii * hbar)) * imag(A[i, j] * H[i, j])
                if gij < 0.0
                    gij = 0.0
                end
                dE = real(H[j, j]) - real(H[i, i])
                if boltz_flag != 0 && dE > eex
                    gij *= exp(-((dE - eex) / (kb * temp)))
                end
                g[i, j] = gij
                g_row_sum += gij
            end
        end
        g[i, i] = 1.0 - g_row_sum
    end
    return nothing
end

@inline function hop_core!(rng::AbstractRNG, g::Matrix{Float64}, hopstate0::Int)
    row = hopstate0 + 1
    ksi = rand(rng)
    nrm = 0.0
    @inbounds for j in 1:size(g, 2)
        nrm += g[row, j]
    end
    if nrm <= 0.0
        return hopstate0
    end
    accum = 0.0
    @inbounds for j in 1:size(g, 2)
        accum += g[row, j] / nrm
        if accum > 0.0 && ksi <= accum
            return j - 1
        end
    end
    return hopstate0
end

@inline function project_out_core!(
    C::Vector{ComplexF64},
    A::Matrix{ComplexF64},
    i0::Int,
    curr_state0::Int,
)
    update_populations_core!(C, A)
    i = i0 + 1
    nrm2 = 0.0
    @inbounds for k in 1:length(C)
        nrm2 += real(A[k, k])
    end
    nrm2 -= real(A[i, i])
    if nrm2 <= 0.0
        zero_state!(C, curr_state0)
    else
        C[i] = 0.0 + 0.0im
        inv = 1.0 / sqrt(nrm2)
        @inbounds for k in 1:length(C)
            C[k] *= inv
        end
    end
    update_populations_core!(C, A)
    return nothing
end

function update_decoherence_times_core!(
    C::Vector{ComplexF64},
    A::Matrix{ComplexF64},
    rates::Matrix{ComplexF64},
    tau_m::Vector{Float64},
)
    update_populations_core!(C, A)
    n = length(C)
    @inbounds for i in 1:n
        s = 0.0
        for j in 1:n
            s += real(rates[i, j]) * real(A[j, j])
        end
        tau_m[i] = s
    end
    return nothing
end

function dish_decoherence_core!(
    rng::AbstractRNG,
    C::Vector{ComplexF64},
    H::AbstractMatrix{ComplexF64},
    A::Matrix{ComplexF64},
    rates::Matrix{ComplexF64},
    tau_m::Vector{Float64},
    t_m::Vector{Float64},
    curr_state0::Int,
    dt::Float64,
    kb::Float64,
    temp::Float64,
)
    update_decoherence_times_core!(C, A, rates, tau_m)
    n = length(C)
    curr = curr_state0
    @inbounds for i0 in 0:(n - 1)
        i = i0 + 1
        if tau_m[i] <= 0.0
            continue
        end
        rnd_i = 1.0 / tau_m[i]
        if t_m[i] >= rnd_i
            zeta = rand(rng)
            P = real(A[i, i])
            dE = real(H[i, i] - H[curr + 1, curr + 1])
            if dE > 0.0
                P *= exp(-(dE / (kb * temp)))
            end
            if zeta < P
                zero_state!(C, i0)
                curr = i0
                break
            else
                project_out_core!(C, A, i0, curr)
                t_m[i] = 0.0
                tau_m[i] = 0.0
            end
        end
    end
    @inbounds for i in 1:n
        t_m[i] += dt
    end
    return curr
end

function sdm_decoherence_core!(
    C::Vector{ComplexF64},
    A::Matrix{ComplexF64},
    rates::Matrix{ComplexF64},
    act_st0::Int,
    dt::Float64,
)
    n = length(C)
    if act_st0 < 0 || act_st0 >= n
        error("Error in sdm_decoherence_core!: active state index out of range")
    end
    act = act_st0 + 1
    update_populations_core!(C, A)
    p_aa_old = real(A[act, act])
    if p_aa_old > 1.0
        sclf = 1.0 / sqrt(p_aa_old)
        @inbounds for i in 1:n
            C[i] *= sclf
        end
        update_populations_core!(C, A)
        p_aa_old = real(A[act, act])
    end
    if p_aa_old <= 0.0
        return nothing
    end
    inact_st_pop = 0.0
    @inbounds for i in 1:n
        if i != act
            sclf = exp(-dt * real(rates[i, act]))
            C[i] *= sclf
            inact_st_pop += abs2(C[i])
        end
    end
    if inact_st_pop > 1.0
        inact_st_pop = 1.0
    end
    p_aa_new = max(0.0, 1.0 - inact_st_pop)
    C[act] *= sqrt(p_aa_new / p_aa_old)
    update_populations_core!(C, A)
    return nothing
end

function gfsh_decohere_core!(
    rng::AbstractRNG,
    C::Vector{ComplexF64},
    H::AbstractMatrix{ComplexF64},
    A::Matrix{ComplexF64},
    tau_m::Vector{Float64},
    t_m::Vector{Float64},
    curr_state0::Int,
    dt::Float64,
    kb::Float64,
    temp::Float64,
)
    update_populations_core!(C, A)
    n = length(C)
    total_pop = 0.0
    @inbounds for i in 1:n
        total_pop += real(A[i, i])
    end
    if total_pop <= 0.0
        return curr_state0
    end
    num_decohere = 0
    @inbounds for i in 1:n
        if real(A[i, i]) / total_pop > rand(rng)
            num_decohere += 1
        end
    end
    if num_decohere <= 0
        @inbounds for i in 1:n
            t_m[i] += dt
        end
        return curr_state0
    end

    perm = collect(0:(n - 1))
    @inbounds for i in n:-1:2
        j = rand(rng, 1:i)
        perm[i], perm[j] = perm[j], perm[i]
    end

    curr = curr_state0
    @inbounds for kk in 1:min(num_decohere, n)
        i0 = perm[kk]
        i = i0 + 1
        zeta = rand(rng)
        P = real(A[i, i]) / total_pop
        dE = real(H[i, i] - H[curr + 1, curr + 1])
        if dE > 0.0
            P *= exp(-(dE / (kb * temp)))
        end
        if zeta < P
            zero_state!(C, i0)
            curr = i0
        else
            project_out_core!(C, A, i0, curr)
        end
        t_m[i] = 0.0
        tau_m[i] = 0.0
    end
    @inbounds for i in 1:n
        t_m[i] += dt
    end
    return curr
end

function build_spin_orbital_H!(
    Htmp::Matrix{ComplexF64},
    Hij::AbstractMatrix{ComplexF64},
    alp_bet::Int,
)
    numstates = size(Hij, 1)
    fill!(Htmp, 0.0 + 0.0im)
    @inbounds for k1 in 0:(numstates - 1), k2 in 0:(numstates - 1)
        i1 = 2 * k1 + 1
        i2 = 2 * k1 + 2
        jx1 = 2 * k2 + 1
        jx2 = 2 * k2 + 2
        val = Hij[k1 + 1, k2 + 1]
        Htmp[i1, jx1] = val
        Htmp[i2, jx2] = val
        if alp_bet != 0
            Htmp[i2, jx1] = val
            Htmp[i1, jx2] = val
        end
    end
    return Htmp
end

function build_Hme_batch(
    H_array::Array{ComplexF64, 3},
    istart::Int,
    namdtime::Int,
    Ham_size::Int,
    me_numstates::Int,
    transition_map::Vector{NTuple{4, Int}},
    diag_orbital_map_packed::Matrix{Int},
    shifts::Vector{Float64},
    alp_bet::Int,
)
    numstates = size(H_array, 2)
    Htmp = zeros(ComplexF64, 2 * numstates, 2 * numstates)
    Hme_batch = zeros(ComplexF64, namdtime, me_numstates, me_numstates)

    @inbounds for j in istart:(istart + namdtime - 1)
        t = j - istart + 1
        j1 = mod(j, Ham_size) + 1
        Hij = @view H_array[j1, :, :]
        build_spin_orbital_H!(Htmp, Hij, alp_bet)

        for (I, J, ii, jj) in transition_map
            Hme_batch[t, I, J] += Htmp[ii, jj]
        end
        for Istate in 1:me_numstates
            diagE = 0.0 + 0.0im
            for k in 1:size(diag_orbital_map_packed, 2)
                orb = diag_orbital_map_packed[Istate, k]
                if orb > 0
                    diagE += Htmp[orb, orb]
                end
            end
            Hme_batch[t, Istate, Istate] = diagE + shifts[Istate]
        end
    end
    return Hme_batch
end

function write_me_energies(
    Hme_batch::Array{ComplexF64, 3},
    istart::Int,
    icond::Int,
    outdir::AbstractString,
)
    if !SAVE_ME_ENERGIES
        return nothing
    end
    outfile = joinpath(outdir, "me_energies$(icond)")
    sz = size(Hme_batch, 1)
    nst = size(Hme_batch, 2)
    open(outfile, "w") do io
        @inbounds for local_t in 1:sz
            j = istart + local_t - 1
            @printf(io, "t= %d  E[0]= %.12g  ", j, real(Hme_batch[local_t, 1, 1]))
            e0 = Hme_batch[local_t, 1, 1]
            for Istate in 0:(nst - 1)
                @printf(
                    io,
                    "E[%d]-E[0]= %.12g  ",
                    Istate,
                    real(Hme_batch[local_t, Istate + 1, Istate + 1] - e0)
                )
            end
            println(io)
        end
    end
end

function autocorr_fft_cross_strict(x::Vector{Float64})
    len = length(x)
    sz = Int(floor(len / 2))
    if sz <= 0
        return Float64[]
    end
    len < 2 * sz &&
        error("Input length is inconsistent with strict autocorrelation definition.")
    a = @view x[1:sz]
    b = @view x[1:(2 * sz - 1)]
    nconv = length(a) + length(b) - 1
    avec = zeros(Float64, nconv)
    bvec = zeros(Float64, nconv)
    @inbounds for i in 1:sz
        avec[i] = a[sz - i + 1]
    end
    @inbounds for i in 1:length(b)
        bvec[i] = b[i]
    end
    conv_full = irfft(rfft(avec) .* rfft(bvec), nconv)
    C = zeros(Float64, sz)
    start_idx = sz
    @inbounds for t in 0:(sz - 1)
        C[t + 1] = conv_full[start_idx + t] / sz
    end
    return C
end

function decoherence_rates(
    x::Vector{Float64},
    dt::Real,
    icond::Int,
    i1::Int,
    j1::Int,
    outdir::AbstractString,
    p::NAMDParams,
)
    len = length(x)
    sz = Int(floor(len / 2))
    C = autocorr_fft_cross_strict(x)
    IC = zeros(Float64, sz)
    IIC = zeros(Float64, sz)

    sum0 = 0.0
    @inbounds for t in 1:sz
        IC[t] = sum0
        sum0 += C[t] * (dt / p.hbar)
    end
    sum0 = 0.0
    @inbounds for t in 1:sz
        IIC[t] = sum0
        sum0 += IC[t] * (dt / p.hbar)
    end

    nrm = C[1]
    Cnorm = copy(C)
    if abs(nrm) > eps(Float64)
        Cnorm ./= nrm
    end

    T = zeros(Float64, sz)
    selIIC = zeros(Float64, sz)
    first_region = true
    cnt = 0
    @inbounds for t in 0:(sz - 1)
        if first_region
            if IIC[t + 1] < 2.3
                cnt += 1
                T[cnt] = t * t * dt * dt
                selIIC[cnt] = IIC[t + 1]
            else
                first_region = false
            end
        end
    end

    a = 0.0
    b = 0.0
    if cnt > 0
        _, _, a, b = regression(@view(T[1:cnt]), @view(selIIC[1:cnt]), a, b)
    end
    if b < 0.0
        b = 0.0
    end

    D = exp.(-IIC)
    dE = 0.0025
    Npoints = 2000
    tt = collect(1:(sz - 1)) .* dt
    Ct = @view Cnorm[2:end]

    spectral_buf = IOBuffer()
    for w in 0:(Npoints - 1)
        ww = w * dE
        s = 0.0
        @inbounds for k in eachindex(Ct)
            s += cos(ww * tt[k]) * Ct[k]
        end
        J = ((1.0 + 2.0 * s) * dt)^2 / (2.0 * pi)
        @printf(
            spectral_buf,
            "w(eV)= %.12g w(cm^-1)= %.12g J= %.12g sqrt(J)= %.12g\n",
            w * dE,
            w * dE * 8065.54468111324,
            J,
            sqrt(J),
        )
    end
    if SAVE_ICOND_FILES
        spectral_path = joinpath(outdir, "icond$(icond)pair$(i1)_$(j1)Spectral_density.txt")
        write(spectral_path, String(take!(spectral_buf)))
    end

    dephase_buf = IOBuffer()
    println(
        dephase_buf,
        "Time    D(t)       fitted D(t)     Normalized_autocorrelation_function  Unnormalized_autocorrelation_function   Second cumulant",
    )
    @inbounds for t in 0:(sz - 1)
        @printf(
            dephase_buf,
            "%.12g  %.12g  %.12g  %.12g %.12g  %.12g\n",
            t * dt,
            D[t + 1],
            exp(-a) * exp(-b * t * t * dt * dt),
            Cnorm[t + 1],
            nrm * Cnorm[t + 1],
            IIC[t + 1],
        )
    end
    if SAVE_ICOND_FILES
        dephase_path =
            joinpath(outdir, "icond$(icond)pair$(i1)_$(j1)Dephasing_function.txt")
        write(dephase_path, String(take!(dephase_buf)))
    end

    return sqrt(b), x
end

function run_decoherence_rates!(
    Hme_batch::Array{ComplexF64, 3},
    icond::Int,
    outdir::AbstractString,
    p::NAMDParams,
)
    sz = size(Hme_batch, 1)
    N = size(Hme_batch, 2)
    rij = zeros(Float64, N, N)

    @threads for pair_idx in 1:(N * N)
        i = div(pair_idx - 1, N) + 1
        j = mod(pair_idx - 1, N) + 1
        if i == j
            rij[i, j] = 0.0
        else
            Eij = zeros(Float64, sz)
            ave_dEij = 0.0
            @inbounds for t in 1:sz
                dEij = real(Hme_batch[t, i, i] - Hme_batch[t, j, j])
                Eij[t] = dEij
                ave_dEij += dEij
            end
            Eij .-= ave_dEij / sz
            rij[i, j], _ = decoherence_rates(Eij, p.nucl_dt, icond, i - 1, j - 1, outdir, p)
        end
    end

    deco_buf = IOBuffer()
    @inbounds for i in 1:N
        for j in 1:N
            @printf(deco_buf, "%.12g ", rij[i, j])
        end
        println(deco_buf)
    end
    if SAVE_DECOHERENCE_RATES
        deco_path = joinpath(outdir, "decoherence_rates_icond$(icond).txt")
        write(deco_path, String(take!(deco_buf)))
    end
    return rij
end

function run_namd_core_single!(
    sh_local::Matrix{Float64},
    se_local::Matrix{Float64},
    Hme_batch::Array{ComplexF64, 3},
    rates::Matrix{ComplexF64},
    init_state::Int,
    traj_start::Int,
    traj_end::Int,
    p::NAMDParams,
    seed_offset::Int,
)
    sz = size(Hme_batch, 1)
    nst = size(Hme_batch, 2)
    nel = Int(round(p.nucl_dt / p.elec_dt))
    rng = MersenneTwister(RANDOM_SEED + seed_offset)
    C = zeros(ComplexF64, nst)
    A = zeros(ComplexF64, nst, nst)
    g = zeros(Float64, nst, nst)
    tau_m = zeros(Float64, nst)
    t_m = zeros(Float64, nst)

    for _ in traj_start:traj_end
        zero_state!(C, init_state)
        curr_state = init_state
        fill!(tau_m, 0.0)
        fill!(t_m, 0.0)
        @inbounds for istep in 1:sz
            H = @view Hme_batch[istep, :, :]
            init_hop_prob_core!(g)
            if p.integrator == 0
                for _el in 1:nel
                    propagate_coefficients_core!(C, H, p.elec_dt, p.hbar)
                    if p.sh_algo == 0
                        update_hop_prob_fssh_core!(
                            C,
                            H,
                            A,
                            g,
                            p.elec_dt,
                            p.hbar,
                            p.kb,
                            p.boltz_flag,
                            p.Temp,
                            0.0,
                        )
                    end
                end
            end
            update_populations_core!(C, A)
            if p.decoherence_method == 0
                curr_state = hop_core!(rng, g, curr_state)
            elseif p.decoherence_method == 1
                curr_state = dish_decoherence_core!(
                    rng,
                    C,
                    H,
                    A,
                    rates,
                    tau_m,
                    t_m,
                    curr_state,
                    p.nucl_dt,
                    p.kb,
                    p.Temp,
                )
            elseif p.decoherence_method == 2
                curr_state = hop_core!(rng, g, curr_state)
                sdm_decoherence_core!(C, A, rates, curr_state, p.nucl_dt)
                update_populations_core!(C, A)
            elseif p.decoherence_method == 3
                curr_state = gfsh_decohere_core!(
                    rng,
                    C,
                    H,
                    A,
                    tau_m,
                    t_m,
                    curr_state,
                    p.nucl_dt,
                    p.kb,
                    p.Temp,
                )
            end
            sh_local[istep, curr_state + 1] += 1.0
            for j in 1:nst
                se_local[istep, j] += real(A[j, j])
            end
        end
    end
    return nothing
end

function run_namd!(
    Hme_batch::Array{ComplexF64, 3},
    icond::Int,
    init_state::Int,
    outdir::AbstractString,
    p::NAMDParams,
    rij_input::Union{Nothing, Matrix{Float64}} = nothing,
)
    sz = size(Hme_batch, 1)
    nst = size(Hme_batch, 2)
    rates = zeros(ComplexF64, nst, nst)
    if p.decoherence_method > 0
        if rij_input !== nothing
            rates .= complex.(rij_input, zeros(Float64, nst, nst))
        else
            filename = joinpath(outdir, "decoherence_rates_icond$(icond).txt")
            r_ij = load_matrix(filename)
            rates .= complex.(r_ij, zeros(Float64, nst, nst))
        end
    end

    nt = min(nthreads(), p.num_sh_traj)
    sh_tls = [zeros(Float64, sz, nst) for _ in 1:nt]
    se_tls = [zeros(Float64, sz, nst) for _ in 1:nt]
    chunks = Vector{Tuple{Int, Int}}(undef, nt)
    for tid in 1:nt
        lo = fld((tid - 1) * p.num_sh_traj, nt) + 1
        hi = fld(tid * p.num_sh_traj, nt)
        chunks[tid] = (lo, hi)
    end

    @threads for tid in 1:nt
        lo, hi = chunks[tid]
        if lo <= hi
            run_namd_core_single!(
                sh_tls[tid],
                se_tls[tid],
                Hme_batch,
                rates,
                init_state,
                lo,
                hi,
                p,
                100000 * icond + tid,
            )
        end
    end

    sh_pops = zeros(Float64, sz, nst)
    se_pops = zeros(Float64, sz, nst)
    for tid in 1:nt
        sh_pops .+= sh_tls[tid]
        se_pops .+= se_tls[tid]
    end
    sh_pops ./= p.num_sh_traj
    se_pops ./= p.num_sh_traj

    out1 = SAVE_ME_POP ? open(joinpath(outdir, "me_pop$(icond)"), "w") : nothing
    out2 = open(joinpath(outdir, "out$(icond)"), "w")

    try
        @inbounds for i in 1:sz
            if out1 !== nothing
                @printf(out1, "time %d ", i - 1)
                tot = 0.0
                for j in 1:nst
                    tot += se_pops[i, j]
                end
                for j in 0:(nst - 1)
                    @printf(out1, "P(%d)= %.10f ", j, se_pops[i, j + 1])
                end
                @printf(out1, "Total= %.10f\n", tot)
            end
            if out2 !== nothing
                @printf(out2, "time %d ", i - 1)
                for j in 0:(nst - 1)
                    @printf(out2, "P(%d)= %.10f ", j, sh_pops[i, j + 1])
                end
                println(out2)
            end
        end
    finally
        if out1 !== nothing
            close(out1)
        end
        close(out2)
    end
    return nothing
end

function load_H_array(rt::AbstractString, Ham_size::Int, active_space::Vector{Int})
    nactive = length(active_space)
    H_array = zeros(ComplexF64, Ham_size, nactive, nactive)
    @inbounds for j in 0:(Ham_size - 1)
        Ham_re = load_matrix(joinpath(rt, "res", "0_Ham_$(j)_re"))
        Ham_im = load_matrix(joinpath(rt, "res", "0_Ham_$(j)_im"))
        Ham =
            complex.(
                extract_2D(Ham_re, active_space, -1),
                extract_2D(Ham_im, active_space, -1),
            ) .* Ry_to_eV
        H_array[j + 1, :, :] .= Ham
    end
    return H_array
end

function main()
    t0 = time()
    filepath = abspath(dirname(PROGRAM_FILE == "" ? (@__FILE__) : PROGRAM_FILE))
    cd(filepath)
    rt = joinpath(filepath, "..")
    outdir = joinpath(filepath, "out")
    isdir(outdir) || mkpath(outdir)

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

    me_numstates = length(states)
    me_states = [make_me_state(s, active_space) for s in states]
    active_map = build_active_map(active_space)
    transition_map = build_transition_map(me_states, active_map)
    diag_orbital_map = build_diag_orbital_map(me_states, active_map)
    diag_orbital_map_packed = pack_diag_orbital_map(diag_orbital_map)
    shifts = [st.shift for st in me_states]

    iconds = zeros(Int, iconds_i, 2)
    for ii in 1:iconds_i
        iconds[ii, :] .= [ii - 1, me_numstates - 1]
    end

    H_array = load_H_array(rt, Ham_size, active_space)

    for icond in 0:(size(iconds, 1) - 1)
        istart = iconds[icond + 1, 1]
        init_state = iconds[icond + 1, 2]
        Hme_batch = build_Hme_batch(
            H_array,
            istart,
            namdtime,
            Ham_size,
            me_numstates,
            transition_map,
            diag_orbital_map_packed,
            shifts,
            alp_bet,
        )
        write_me_energies(Hme_batch, istart, icond, outdir)
        rij = nothing
        if p.decoherence_method > 0
            rij = run_decoherence_rates!(Hme_batch, icond, outdir, p)
        end
        run_namd!(Hme_batch, icond, init_state, outdir, p, rij)
        Hme_batch = nothing
        GC.gc()
    end

    elapsedTime = time() - t0
    open(joinpath(filepath, "log.txt"), "w") do logFile
        @printf("Running time: %.1f seconds\n", elapsedTime)
        @printf(logFile, "Running time: %.1f seconds\n", elapsedTime)
    end

    mepop0 = joinpath(outdir, "me_pop0")
    if SAVE_ME_POP && isfile(mepop0)
        last_line = ""
        open(mepop0, "r") do io
            for line in eachline(io)
                last_line = line
            end
        end
        println(last_line)
    end
end

if PROGRAM_FILE == "" || abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
