using LinearAlgebra
using DelimitedFiles
using Printf
using Random
using Statistics

const RANDOM_SEED = 42
Random.seed!(RANDOM_SEED)

global active_space = Int[]
global active_map = Dict{Int, Int}()
global hbar = 0.658218
global kb = 8.617e-5
global Ry_to_eV = 13.60569253
global Ha_to_eV = 27.211396
global rt = ""
global outdir = ""
global nucl_dt = 1.0
global elec_dt = 1.0
global integrator = 0
global sh_algo = 0
global num_sh_traj = 1000
global decoherence_method = 2
global boltz_flag = 1
global Temp = 300.0

mutable struct ElectronicStructure
    curr_state::Int
    num_states::Int
    Ccurr::Vector{ComplexF64}
    Cprev::Vector{ComplexF64}
    Cnext::Vector{ComplexF64}
    Hcurr::Matrix{ComplexF64}
    Hprev::Matrix{ComplexF64}
    Hnext::Matrix{ComplexF64}
    dHdt::Matrix{ComplexF64}
    g::Matrix{Float64}
    tau_m::Vector{Float64}
    t_m::Vector{Float64}
    A::Matrix{ComplexF64}
end

function ElectronicStructure(n::Integer)
    n = Int(n)
    return ElectronicStructure(
        0,
        n,
        zeros(ComplexF64, n),
        zeros(ComplexF64, n),
        zeros(ComplexF64, n),
        zeros(ComplexF64, n, n),
        zeros(ComplexF64, n, n),
        zeros(ComplexF64, n, n),
        zeros(ComplexF64, n, n),
        Matrix{Float64}(I, n, n),
        zeros(Float64, n),
        zeros(Float64, n),
        zeros(ComplexF64, n, n),
    )
end

mutable struct MEState
    name::String
    active_space::Vector{Int}
    actual_state::Vector{Int}
    shift::Float64
    nac_scl_indx::Vector{Int}
    nac_scl::Vector{Float64}
end

function make_me_state(cs)
    global active_space
    shift = length(cs) >= 3 ? Float64(cs[3]) : 0.0
    return MEState(cs[1], copy(active_space), copy(cs[2]), shift, Int[], Float64[])
end

function load_matrix(path::AbstractString)
    if !isfile(path)
        error("Cannot find file: $path")
    end
    data = readdlm(path)
    return Matrix{Float64}(data)
end

function extract_2D(input_matrix::AbstractMatrix, templ::Vector{Int}, shift::Int)
    idx = templ .+ shift .+ 1
    return input_matrix[idx, idx]
end

function ext2int(external::Int, active_space_arg::Vector{Int})
    global active_map
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
    if abs(sx) < eps(Float64)
        b = 0.0
    else
        b = sum(Y) / sx
    end
    return X, Y, a, b
end

function set_state!(es::ElectronicStructure, indx::Int)
    es.Ccurr .= 0.0 + 0.0im
    es.Ccurr[indx + 1] = 1.0 + 0.0im
    es.curr_state = indx
    return es
end

function insertion!(this::ElectronicStructure, es::ElectronicStructure)
    this.num_states = es.num_states
    this.curr_state = es.curr_state
    this.Cprev .= es.Cprev
    this.Ccurr .= es.Ccurr
    this.Cnext .= es.Cnext
    this.tau_m .= es.tau_m
    this.t_m .= es.t_m
    this.A .= es.A
    return this
end

function init_hop_prob1!(es::ElectronicStructure)
    es.g .= Matrix{Float64}(I, es.num_states, es.num_states)
    return es
end

function Efield(t::Real, E::Vector{Float64}, Eex::Real)
    E .= 0.0
    return E, 0.0
end

function rot1!(es::ElectronicStructure, phi::Real, i::Int, j::Int)
    c = cos(phi)
    s = sin(phi)
    ci = es.Ccurr[i + 1]
    cj = es.Ccurr[j + 1]
    es.Ccurr[i + 1] = c * ci + s * cj
    es.Ccurr[j + 1] = -s * ci + c * cj
    return es
end

function rot2!(es::ElectronicStructure, phi::Real, i::Int, j::Int)
    cs = cos(phi)
    isi = 1im * sin(phi)
    ci = es.Ccurr[i + 1]
    cj = es.Ccurr[j + 1]
    es.Ccurr[i + 1] = cs * ci + isi * cj
    es.Ccurr[j + 1] = isi * ci + cs * cj
    return es
end

function rot!(es::ElectronicStructure, Hij::Complex, dt::Real, i::Int, j::Int)
    global hbar
    phi1 = 0.5 * dt * imag(Hij) / hbar
    phi2 = -dt * real(Hij) / hbar
    rot1!(es, phi1, i, j)
    rot2!(es, phi2, i, j)
    rot1!(es, phi1, i, j)
    return es
end

function phase!(es::ElectronicStructure, Hii::Complex, dt::Real, i::Int)
    global hbar
    phi = -dt * real(Hii) / hbar
    es.Ccurr[i + 1] = exp(1im * phi) * es.Ccurr[i + 1]
    return es
end

function propagate_coefficients!(es::ElectronicStructure, dt::Real, Ef::Vector{Float64})
    for i in 0:(es.num_states - 1)
        for j in (i + 1):(es.num_states - 1)
            rot!(es, es.Hcurr[i + 1, j + 1], 0.5 * dt, i, j)
        end
    end
    for i in 0:(es.num_states - 1)
        phase!(es, es.Hcurr[i + 1, i + 1], dt, i)
    end
    for i in (es.num_states - 1):-1:0
        for j in (es.num_states - 1):-1:(i + 1)
            rot!(es, es.Hcurr[i + 1, j + 1], 0.5 * dt, i, j)
        end
    end
    return es
end

function update_populations!(es::ElectronicStructure)
    es.A .= conj.(es.Ccurr) * transpose(es.Ccurr)
    return es
end

function update_hop_prob_fssh!(
    es::ElectronicStructure,
    dt::Real,
    boltz_flag_arg::Int,
    Temp_arg::Real,
    Ef::Vector{Float64},
    Eex::Real,
    rates::Matrix{ComplexF64},
)
    global hbar, kb
    update_populations!(es)
    Heff = es.Hcurr
    for i in 0:(es.num_states - 1)
        a_ii = real(es.A[i + 1, i + 1])
        if a_ii < 1e-12
            a_ii = 1e-12
        end
        g_row_sum = 0.0
        for j in 0:(es.num_states - 1)
            if j != i
                g_ij =
                    (2.0 * dt / (a_ii * hbar)) *
                    imag(es.A[i + 1, j + 1] * Heff[i + 1, j + 1])
                if g_ij < 0.0
                    g_ij = 0.0
                end
                E_i = real(Heff[i + 1, i + 1])
                E_j = real(Heff[j + 1, j + 1])
                dE = E_j - E_i
                bf = 1.0
                if dE > Eex
                    bf = exp(-((dE - Eex) / (kb * Temp_arg)))
                end
                g_ij *= bf
                g_row_sum += g_ij
                es.g[i + 1, j + 1] = g_ij
            end
        end
        es.g[i + 1, i + 1] = 1.0 - g_row_sum
    end
    return es
end

function propagate_electronic!(
    es::Vector{ElectronicStructure},
    i::Int,
    rates::Matrix{ComplexF64},
)
    global integrator, sh_algo, nucl_dt, elec_dt, boltz_flag, Temp
    nel = Int(round(nucl_dt / elec_dt))
    Eex = 0.0
    Ef = zeros(Float64, 3)
    if integrator == 0
        for j in 0:(nel - 1)
            tim = i * nucl_dt + j * elec_dt
            Ef, Eex = Efield(tim, Ef, Eex)
            propagate_coefficients!(es[i + 1], elec_dt, Ef)
            if sh_algo == 0
                update_hop_prob_fssh!(es[i + 1], elec_dt, boltz_flag, Temp, Ef, Eex, rates)
            end
        end
    end
    return es, rates
end

function hop(sh_prob::Matrix{Float64}, hopstate::Int, numstates::Int)
    input_state = hopstate
    ksi = rand()
    probs = vec(sh_prob[input_state + 1, :])
    nrm = sum(probs)
    if nrm <= 0
        return sh_prob, hopstate
    end
    probs ./= nrm
    cprobs = cumsum(probs)
    idx = findfirst(x -> (0 < x && ksi <= x), cprobs)
    hopstate = isnothing(idx) ? input_state : idx - 1
    return sh_prob, hopstate
end

function decoherence_rates(x::Vector{Float64}, dt::Real, icond::Int, i1::Int, j1::Int)
    global hbar, outdir
    len = length(x)
    sz = Int(floor(len / 2))
    C = zeros(Float64, sz)
    IC = zeros(Float64, sz)
    IIC = zeros(Float64, sz)
    D = zeros(Float64, sz)
    T = zeros(Float64, sz)
    selIIC = zeros(Float64, sz)

    for t in 0:(sz - 1)
        C[t + 1] = sum(x[1:sz] .* x[(1 + t):(sz + t)]) / sz
    end

    sum0 = 0.0
    for t in 0:(sz - 1)
        IC[t + 1] = sum0
        sum0 += C[t + 1] * (dt / hbar)
    end

    sum0 = 0.0
    for t in 0:(sz - 1)
        IIC[t + 1] = sum0
        sum0 += IC[t + 1] * (dt / hbar)
    end

    D .= exp.(-IIC)
    nrm = C[1]
    if abs(nrm) > eps(Float64)
        C ./= nrm
    end

    dE = 0.0025
    Npoints = 2000
    J = zeros(Float64, Npoints)
    tt = collect(1:(sz - 1)) .* dt
    Ct = C[2:end]

    spectral_path = joinpath(outdir, "icond$(icond)pair$(i1)_$(j1)Spectral_density.txt")
    open(spectral_path, "w") do io
        for w in 0:(Npoints - 1)
            ww = w * dE
            J[w + 1] = (1.0 + 2.0 * sum(cos.(ww .* tt) .* Ct)) * dt
            J[w + 1] = J[w + 1]^2 / (2.0 * pi)
            @printf(
                io,
                "w(eV)= %.12g w(cm^-1)= %.12g J= %.12g sqrt(J)= %.12g\n",
                w * dE,
                w * dE * 8065.54468111324,
                J[w + 1],
                sqrt(J[w + 1])
            )
        end
    end

    first_region = true
    cnt = 0
    for t in 0:(sz - 1)
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

    T = T[1:cnt]
    selIIC = selIIC[1:cnt]
    a = 0.0
    b = 0.0
    if !isempty(T)
        _, _, a, b = regression(T, selIIC, a, b)
    end
    if b < 0.0
        b = 0.0
    end

    dephase_path = joinpath(outdir, "icond$(icond)pair$(i1)_$(j1)Dephasing_function.txt")
    open(dephase_path, "w") do io
        println(
            io,
            "Time    D(t)       fitted D(t)     Normalized_autocorrelation_function  Unnormalized_autocorrelation_function   Second cumulant",
        )
        for t in 0:(sz - 1)
            @printf(
                io,
                "%.12g  %.12g  %.12g  %.12g %.12g  %.12g\n",
                t * dt,
                D[t + 1],
                exp(-a) * exp(-b * t * t * dt * dt),
                C[t + 1],
                nrm * C[t + 1],
                IIC[t + 1]
            )
        end
    end

    return sqrt(b), x
end

function run_decoherence_rates!(me_es::Vector{ElectronicStructure}, icond::Int)
    global nucl_dt, outdir
    sz = length(me_es)
    N = me_es[1].num_states
    rij = zeros(Float64, N, N)
    deco_path = joinpath(outdir, "decoherence_rates_icond$(icond).txt")
    open(deco_path, "w") do io
        for i1 in 0:(N - 1)
            for j1 in 0:(N - 1)
                if i1 == j1
                    rij[i1 + 1, j1 + 1] = 0.0
                else
                    Eij = zeros(Float64, sz)
                    ave_dEij = 0.0
                    for t in 0:(sz - 1)
                        dEij = real(
                            me_es[t + 1].Hcurr[i1 + 1, i1 + 1] -
                            me_es[t + 1].Hcurr[j1 + 1, j1 + 1],
                        )
                        Eij[t + 1] = dEij
                        ave_dEij += dEij
                    end
                    Eij .-= ave_dEij / sz
                    rij[i1 + 1, j1 + 1], _ = decoherence_rates(Eij, nucl_dt, icond, i1, j1)
                end
                @printf(io, "%.12g ", real(rij[i1 + 1, j1 + 1]))
            end
            println(io)
        end
    end
    return me_es
end

function update_decoherence_times!(es::ElectronicStructure, rates::Matrix{ComplexF64})
    update_populations!(es)
    pops = real.(diag(es.A))
    es.tau_m .= real.(rates) * pops
    return es
end

function decohere!(es::ElectronicStructure, i::Int)
    es.Ccurr .= 0.0 + 0.0im
    es.Ccurr[i + 1] = 1.0 + 0.0im
    es.curr_state = i
    return es
end

function project_out!(es::ElectronicStructure, i::Int)
    update_populations!(es)
    es.Ccurr[i + 1] = 0.0 + 0.0im
    nrm2 = sum(real.(diag(es.A))) - real(es.A[i + 1, i + 1])
    if nrm2 <= 0
        es.Ccurr .= 0.0 + 0.0im
        es.Ccurr[es.curr_state + 1] = 1.0 + 0.0im
    else
        es.Ccurr ./= sqrt(nrm2)
    end
    update_populations!(es)
    return es
end

function dish_decoherence!(
    es::ElectronicStructure,
    dt::Real,
    boltz_flag_arg::Int,
    Temp_arg::Real,
    rates::Matrix{ComplexF64},
)
    global kb
    update_decoherence_times!(es, rates)
    for i in 0:(es.num_states - 1)
        if es.tau_m[i + 1] <= 0
            continue
        end
        rnd_i = 1.0 / es.tau_m[i + 1]
        if es.t_m[i + 1] >= rnd_i
            zeta = rand()
            P = real(es.A[i + 1, i + 1])
            dE = real(
                es.Hcurr[i + 1, i + 1] - es.Hcurr[es.curr_state + 1, es.curr_state + 1],
            )
            if dE > 0
                P *= exp(-(dE / (kb * Temp_arg)))
            end
            if zeta < P
                decohere!(es, i)
                break
            else
                project_out!(es, i)
                es.t_m[i + 1] = 0.0
                es.tau_m[i + 1] = 0.0
            end
        end
    end
    es.t_m .+= dt
    return es
end

function sdm_decoherence!(
    es::ElectronicStructure,
    dt::Real,
    act_st::Int,
    rates::Matrix{ComplexF64},
    tol::Real = 0.0,
)
    if act_st < 0 || act_st >= es.num_states
        error(
            "Error in ElectronicStructure_sdm_decoherence: active state index out of range",
        )
    end
    update_populations!(es)
    p_aa_old = real(es.A[act_st + 1, act_st + 1])
    if p_aa_old > 1.0 + tol
        sclf = 1.0 / sqrt(p_aa_old)
        es.Ccurr .*= sclf
        update_populations!(es)
        p_aa_old = real(es.A[act_st + 1, act_st + 1])
    end
    if p_aa_old <= 0.0
        return es
    end
    inact_st_pop = 0.0
    for i in 0:(es.num_states - 1)
        if i != act_st
            itau = real(rates[i + 1, act_st + 1])
            sclf = exp(-dt * itau)
            es.Ccurr[i + 1] *= sclf
            inact_st_pop += abs(es.Ccurr[i + 1])^2
        end
    end
    if inact_st_pop > 1.0
        error(
            "Error in ElectronicStructure_sdm_decoherence: inactive-state population > 1.0",
        )
    end
    p_aa_new = 1.0 - inact_st_pop
    if p_aa_new < 0.0
        error(
            "Error in ElectronicStructure_sdm_decoherence: new active-state population is negative",
        )
    end
    sclf = sqrt(p_aa_new / p_aa_old)
    es.Ccurr[act_st + 1] *= sclf
    update_populations!(es)
    new_norm = real(dot(es.Ccurr, es.Ccurr))
    if abs(new_norm - 1.0) > 0.1
        error("Error in ElectronicStructure_sdm_decoherence: norm deviates too much from 1")
    end
    return es
end

function gfsh_decohere!(es::ElectronicStructure, dt::Real, Temp_arg::Real)
    global kb
    update_populations!(es)
    curr_state = es.curr_state
    pop_diag = real.(diag(es.A))
    total_pop = sum(pop_diag)
    if total_pop <= 0
        return es
    end
    normalized_pop = pop_diag ./ total_pop
    num_decohere = sum(normalized_pop .> rand(length(normalized_pop)))
    if num_decohere <= 0
        es.t_m .+= dt
        return es
    end
    decohere_states = randperm(es.num_states)[1:min(num_decohere, es.num_states)] .- 1
    for i in decohere_states
        zeta = rand()
        P = real(es.A[i + 1, i + 1]) / total_pop
        dE = real(es.Hcurr[i + 1, i + 1] - es.Hcurr[curr_state + 1, curr_state + 1])
        if dE > 0
            P *= exp(-(dE / (kb * Temp_arg)))
        end
        if zeta < P
            decohere!(es, i)
            curr_state = i
        else
            project_out!(es, i)
        end
        es.t_m[i + 1] = 0.0
        es.tau_m[i + 1] = 0.0
    end
    es.t_m .+= dt
    return es
end

# =========================
# NAMD main routine
# =========================
function run_namd!(
    me_es::Vector{ElectronicStructure},
    me_states::Vector{MEState},
    icond::Int,
)
    global nucl_dt, elec_dt, num_sh_traj, decoherence_method, outdir, Temp
    sz = length(me_es)
    nst = me_es[1].num_states
    init_state = me_es[1].curr_state
    sh_pops = zeros(Float64, sz, nst)
    se_pops = zeros(Float64, sz, nst)
    rates = zeros(ComplexF64, nst, nst)

    if decoherence_method > 0
        filename = joinpath(outdir, "decoherence_rates_icond$(icond).txt")
        r_ij = load_matrix(filename)
        rates .= complex.(r_ij, zeros(Float64, nst, nst))
    end

    for n in 1:num_sh_traj
        set_state!(me_es[1], init_state)
        me_es[1].t_m[1] = 0.0
        for i in 0:(sz - 1)
            if i > 0
                insertion!(me_es[i + 1], me_es[i])
            end
            init_hop_prob1!(me_es[i + 1])
            me_es, rates = propagate_electronic!(me_es, i, rates)
            update_populations!(me_es[i + 1])

            if decoherence_method == 0
                me_es[i + 1].g, me_es[i + 1].curr_state =
                    hop(me_es[i + 1].g, me_es[i + 1].curr_state, nst)
            elseif decoherence_method == 1
                dish_decoherence!(me_es[i + 1], nucl_dt, 1, Temp, rates)
            elseif decoherence_method == 2
                me_es[i + 1].g, me_es[i + 1].curr_state =
                    hop(me_es[i + 1].g, me_es[i + 1].curr_state, nst)
                curr_state = me_es[i + 1].curr_state
                sdm_decoherence!(me_es[i + 1], nucl_dt, curr_state, rates)
                update_populations!(me_es[i + 1])
            elseif decoherence_method == 3
                gfsh_decohere!(me_es[i + 1], nucl_dt, Temp)
            end

            curr_state = me_es[i + 1].curr_state
            sh_pops[i + 1, curr_state + 1] += 1.0
            se_pops[i + 1, :] .+= vec(real.(diag(me_es[i + 1].A)))
        end
    end

    se_pops2 = se_pops ./ num_sh_traj
    sh_pops2 = sh_pops ./ num_sh_traj
    outfile1 = joinpath(outdir, "me_pop$(icond)")
    outfile2 = joinpath(outdir, "out$(icond)")

    open(outfile1, "w") do out1
        open(outfile2, "w") do out2
            for i in 0:(sz - 1)
                @printf(out1, "time %d ", i)
                tot = sum(se_pops2[i + 1, :])
                for j in 0:(nst - 1)
                    @printf(out1, "P(%d)= %.10f ", j, se_pops2[i + 1, j + 1])
                end
                @printf(out1, "Total= %.10f\n", tot)

                @printf(out2, "time %d ", i)
                for j in 0:(nst - 1)
                    @printf(out2, "P(%d)= %.10f ", j, sh_pops2[i + 1, j + 1])
                end
                println(out2)
            end
        end
    end

    return me_es, me_states
end

function main()
    global active_space, hbar, kb, Ry_to_eV, Ha_to_eV, rt, outdir
    global nucl_dt, elec_dt, integrator, sh_algo, num_sh_traj
    global decoherence_method, boltz_flag, Temp, active_map

    t0 = time()

    filepath = abspath(dirname(PROGRAM_FILE == "" ? (@__FILE__) : PROGRAM_FILE))
    cd(filepath)
    rt = joinpath(filepath, "..")
    outdir = joinpath(filepath, "out")
    if !isdir(outdir)
        mkpath(outdir)
    end

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
    hbar = 0.658218
    kb = 8.617e-5
    Ry_to_eV = 13.60569253
    Ha_to_eV = 27.211396

    me_numstates = length(states)
    numstates = length(active_space)
    me_states = [make_me_state(s) for s in states]
    num_elec = length(me_states[1].actual_state)

    active_map = Dict{Int, Int}()
    for k in 1:numstates
        active_map[active_space[k]] = k - 1
    end

    iconds = zeros(Int, iconds_i, 2)
    for ii in 1:iconds_i
        iconds[ii, :] .= [ii - 1, me_numstates - 1]
    end

    H_batch = Vector{Matrix{ComplexF64}}(undef, Ham_size)
    for j in 0:(Ham_size - 1)
        Ham_re = load_matrix(joinpath(rt, "res", "0_Ham_$(j)_re"))
        Ham_im = load_matrix(joinpath(rt, "res", "0_Ham_$(j)_im"))
        Ham = complex.(
            extract_2D(Ham_re, active_space, -1),
            extract_2D(Ham_im, active_space, -1),
        )
        H_batch[j + 1] = Ham .* Ry_to_eV
    end

    for icond in 0:(size(iconds, 1) - 1)
        oe_es = [ElectronicStructure(2 * numstates) for _ in 1:namdtime]
        me_es = [ElectronicStructure(me_numstates) for _ in 1:namdtime]

        istart = iconds[icond + 1, 1]
        for j in istart:(istart + namdtime - 1)
            t = j - istart + 1
            j1 = mod(j, Ham_size)
            Hij = H_batch[j1 + 1]
            Htmp = zeros(ComplexF64, 2 * numstates, 2 * numstates)
            for k1 in 0:(numstates - 1)
                i1 = 2 * k1 + 1
                i2 = 2 * k1 + 2
                for k2 in 0:(numstates - 1)
                    jx1 = 2 * k2 + 1
                    jx2 = 2 * k2 + 2
                    Htmp[i1, jx1] = Hij[k1 + 1, k2 + 1]
                    Htmp[i2, jx2] = Hij[k1 + 1, k2 + 1]
                end
            end

            if alp_bet != 0
                for k1 in 0:(numstates - 1)
                    i1 = 2 * k1 + 1
                    i2 = 2 * k1 + 2
                    for k2 in 0:(numstates - 1)
                        jx1 = 2 * k2 + 1
                        jx2 = 2 * k2 + 2
                        Htmp[i2, jx1] = Hij[k1 + 1, k2 + 1]
                        Htmp[i1, jx2] = Hij[k1 + 1, k2 + 1]
                    end
                end
            end
            oe_es[t].Hcurr .= Htmp
        end

        set_state!(me_es[1], iconds[icond + 1, 2])

        for j in istart:(istart + namdtime - 1)
            t = j - istart + 1
            Hme = zeros(ComplexF64, me_numstates, me_numstates)
            for Istate in 0:(me_numstates - 1)
                for Jstate in 0:(me_numstates - 1)
                    orb_i = 0
                    orb_j = 0
                    delt, _, _, orb_i, orb_j = delta_states(
                        me_states[Istate + 1].actual_state,
                        me_states[Jstate + 1].actual_state,
                        orb_i,
                        orb_j,
                    )
                    if delt
                        orb_i = ext2int(orb_i, active_space)
                        orb_j = ext2int(orb_j, active_space)
                        Hme[Istate + 1, Jstate + 1] += oe_es[t].Hcurr[orb_i + 1, orb_j + 1]
                    end
                end

                diagE = 0.0 + 0.0im
                for el in 0:(num_elec - 1)
                    orb_i =
                        ext2int(me_states[Istate + 1].actual_state[el + 1], active_space)
                    diagE += oe_es[t].Hcurr[orb_i + 1, orb_i + 1]
                end
                shift = me_states[Istate + 1].shift
                Hme[Istate + 1, Istate + 1] = diagE + shift
            end
            me_es[t].Hcurr .= Hme
        end

        outfile = joinpath(outdir, "me_energies$(icond)")
        open(outfile, "w") do io
            for j in istart:(istart + namdtime - 1)
                t = j - istart + 1
                @printf(io, "t= %d  E[0]= %.12g  ", j, real(me_es[t].Hcurr[1, 1]))
                e0 = me_es[t].Hcurr[1, 1]
                for Istate in 0:(me_es[t].num_states - 1)
                    @printf(
                        io,
                        "E[%d]-E[0]= %.12g  ",
                        Istate,
                        real(me_es[t].Hcurr[Istate + 1, Istate + 1] - e0)
                    )
                end
                println(io)
            end
        end

        if decoherence_method > 0
            run_decoherence_rates!(me_es, icond)
        end

        me_es, me_states = run_namd!(me_es, me_states, icond)
        oe_es = nothing
        me_es = nothing
        GC.gc()
    end

    elapsedTime = time() - t0
    open(joinpath(filepath, "log.txt"), "w") do logFile
        @printf("Running time: %.1f seconds\n", elapsedTime)
        @printf(logFile, "Running time: %.1f seconds\n", elapsedTime)
    end

    mepop0 = joinpath(outdir, "me_pop0")
    if isfile(mepop0)
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
