import Pkg
using Dates
using LinearAlgebra
using Printf
using Random
using Statistics

VERSION >= v"1.10.0" || error("Julia 1.10 or newer is required.")

const ROOT = @__DIR__
const ENV_DIR = joinpath(ROOT, ".julia_env")
const REQUIRED_DFTK = v"0.7.26"
const AUTO_INSTALL = lowercase(get(ENV, "ADAPTIVEEXX_AUTO_INSTALL", "1")) in ("1", "true", "yes", "on")
const PSP_NAME = "dojo.nc.sr.pbe.v0_4_1.standard.upf"
const ANG2BOHR = 1.8897261254578281
const HARTREE2EV = 27.211386245988
const MAIN_DIR = joinpath(ROOT, "data", "main")
const SI_DIR = joinpath(ROOT, "data", "si")
const FIG_DIR = joinpath(ROOT, "figures")
const LOG_FILE = joinpath(ROOT, "log.dat")
const RNG_SEED = 20260811
const DFTK_SEED = UInt64(RNG_SEED)

mkpath(MAIN_DIR)
mkpath(SI_DIR)
mkpath(FIG_DIR)

function depinfo(name)
    for dep in values(Pkg.dependencies())
        dep.name == name && return dep
    end
    nothing
end

function depversion(name)
    info = depinfo(name)
    info === nothing ? nothing : info.version
end

function is_direct_dependency(name)
    info = depinfo(name)
    info !== nothing && info.is_direct_dep
end

function ensure_package(name; version=nothing)
    info = depinfo(name)
    current = info === nothing ? nothing : info.version
    direct = info !== nothing && info.is_direct_dep
    if version === nothing
        direct && return
        AUTO_INSTALL || error("$name is not a direct dependency of the active project. Run: import Pkg; Pkg.add(\"$name\")")
        Pkg.add(Pkg.PackageSpec(name=name))
    else
        direct && current == version && return
        AUTO_INSTALL || error("$name $version must be a direct dependency of the active project. Run: import Pkg; Pkg.add(Pkg.PackageSpec(name=\"$name\", version=\"$version\"))")
        if info !== nothing && info.is_pinned
            try
                Pkg.free(name)
            catch
            end
        end
        Pkg.add(Pkg.PackageSpec(name=name, version=version))
    end
end

function environment_ready()
    dftk = depinfo("DFTK")
    psp = depinfo("PseudoPotentialData")
    dftk !== nothing || return false
    psp !== nothing || return false
    dftk.is_direct_dep || return false
    psp.is_direct_dep || return false
    dftk.version == REQUIRED_DFTK || return false
    dftk.is_pinned || return false
    Base.find_package("DFTK") !== nothing || return false
    Base.find_package("PseudoPotentialData") !== nothing || return false
    true
end

Pkg.activate(ENV_DIR; shared=false)
if environment_ready()
    println("Julia environment ready; skipping package resolution and instantiation.")
else
    println("Preparing pinned Julia environment for first use or repair.")
    ensure_package("DFTK"; version=REQUIRED_DFTK)
    Pkg.pin(Pkg.PackageSpec(name="DFTK"))
    ensure_package("PseudoPotentialData")
    Pkg.resolve()
    Pkg.instantiate()
end

using DFTK
using PseudoPotentialData

Base.pkgversion(DFTK) == REQUIRED_DFTK || error("DFTK version mismatch: $(Base.pkgversion(DFTK)); required $REQUIRED_DFTK")

const PSP = PseudoFamily(PSP_NAME)
const PBE_CACHE = Dict{Any,Any}()

function nowstamp()
    Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS")
end

function write_log_header()
    open(LOG_FILE, "w") do io
        println(io, "project\tAdaptiveEXX.jl")
        println(io, "timestamp\t$(nowstamp())")
        println(io, "julia_version\t$(VERSION)")
        println(io, "dftk_version\t$(Base.pkgversion(DFTK))")
        println(io, "pseudopotentialdata_version\t$(depversion("PseudoPotentialData"))")
        println(io, "pseudopotential_family\t$PSP_NAME")
        println(io, "backend\tCPU")
        println(io, "cpu\t$(Sys.CPU_NAME)")
        println(io, "os\t$(Sys.KERNEL)")
        println(io, "architecture\t$(Sys.ARCH)")
        println(io, "cpu_threads\t$(Sys.CPU_THREADS)")
        println(io, "julia_threads\t$(Threads.nthreads())")
        println(io, "blas_threads\t$(BLAS.get_num_threads())")
        println(io, "hybrid_kpoint_scope\tGamma only")
        println(io, "total_memory_GiB\t$(round(Sys.total_memory()/2.0^30, digits=3))")
        println(io, "free_memory_GiB\t$(round(Sys.free_memory()/2.0^30, digits=3))")
        println(io, "rng_seed\t$RNG_SEED")
        println(io, "dftk_seed_type\t$(typeof(DFTK_SEED))")
        println(io, "controller_policy\tsingle-probe stochastic defect monitoring with defect-triggered rebuild and five-call safety cap")
    end
end

function append_log(fields...)
    open(LOG_FILE, "a") do io
        println(io, join(string.(fields), '\t'))
    end
end

function write_tsv(path, header, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, '\t'))
        for row in rows
            println(io, join(string.(row), '\t'))
        end
    end
end

function diamond(sym, a_ang)
    a = a_ang * ANG2BOHR
    lattice = (a / 2) .* [0.0 1.0 1.0; 1.0 0.0 1.0; 1.0 1.0 0.0]
    X = ElementPsp(sym, PSP)
    lattice, [X, X], [ones(3)/8, -ones(3)/8]
end

function rocksalt(a_sym, b_sym, a_ang)
    a = a_ang * ANG2BOHR
    lattice = (a / 2) .* [0.0 1.0 1.0; 1.0 0.0 1.0; 1.0 1.0 0.0]
    A = ElementPsp(a_sym, PSP)
    B = ElementPsp(b_sym, PSP)
    lattice, [A, B], [zeros(3), fill(0.5, 3)]
end

function replicate_structure(lattice, atoms, positions, reps)
    nx, ny, nz = reps
    newlat = lattice * Diagonal([nx, ny, nz])
    newatoms = eltype(atoms)[]
    newpos = Vector{Vector{Float64}}()
    for i in 0:nx-1, j in 0:ny-1, k in 0:nz-1
        shift = [i, j, k]
        for (a, pos) in zip(atoms, positions)
            push!(newatoms, a)
            push!(newpos, (collect(pos) .+ shift) ./ [nx, ny, nz])
        end
    end
    newlat, newatoms, newpos
end

function systems()
    si2 = diamond(:Si, 5.431)
    c2 = diamond(:C, 3.567)
    mgo = rocksalt(:Mg, :O, 4.212)
    lif = rocksalt(:Li, :F, 4.027)
    si8 = replicate_structure(si2..., (2,2,1))
    [
        (name="Si2", geom=si2, Ecut=38.0, functional="PBE0"),
        (name="C2", geom=c2, Ecut=48.0, functional="PBE0"),
        (name="MgO2", geom=mgo, Ecut=48.0, functional="PBE0"),
        (name="LiF2", geom=lif, Ecut=55.0, functional="PBE0"),
        (name="Si8", geom=si8, Ecut=34.0, functional="PBE0"),
    ]
end

function si16_system()
    si2 = diamond(:Si, 5.431)
    (name="Si16", geom=replicate_structure(si2..., (2,2,2)), Ecut=34.0, functional="PBE0")
end

function hse_systems()
    [
        (name="Si2_HSE06", geom=diamond(:Si, 5.431), Ecut=38.0, functional="HSE06"),
        (name="MgO2_HSE06", geom=rocksalt(:Mg, :O, 4.212), Ecut=48.0, functional="HSE06"),
    ]
end

function occupied_exchange_data(basis, kpt, interaction_kernel, ψk, occk, maskk_occ)
    ψk_occ = @view ψk[:, maskk_occ]
    ψk_occ_real = similar(ψk, basis.fft_size..., length(maskk_occ))
    for (ψnk_real, ψnk) in zip(eachslice(ψk_occ_real; dims=4), eachcol(ψk_occ))
        DFTK.ifft!(ψnk_real, basis, kpt, ψnk)
    end
    occk_occ = @view occk[maskk_occ]
    Vxk = DFTK.ExchangeOperator(basis, kpt, interaction_kernel, occk_occ, ψk_occ, ψk_occ_real)
    Vxk, ψk_occ, ψk_occ_real, occk_occ
end

function exact_action(Vxk, ψ)
    Vxk * ψ
end

ace_action(op, ψ) = op.P * (op.D * (op.P' * ψ))

function full_occupied_defect(Vxk, cached_op, ψk, maskk_occ)
    isempty(maskk_occ) && return 0.0
    num = 0.0
    den = 0.0
    for i in maskk_occ
        v = @view ψk[:, i]
        exact = exact_action(Vxk, v)
        approx = ace_action(cached_op, v)
        num += norm(exact - approx)^2
        den += norm(exact)^2
    end
    sqrt(num / max(den, eps()))
end

function stochastic_defect(Vxk, cached_op, ψk, maskk_occ, nprobe, seed)
    nocc = length(maskk_occ)
    (nocc == 0 || nprobe <= 0) && return (0.0, 0.0, 0.0)
    Ψ = @view ψk[:, maskk_occ]
    rng = MersenneTwister(seed)
    vals = Float64[]
    num = 0.0
    den = 0.0
    for _ in 1:nprobe
        z = [rand(rng, Bool) ? 1.0 : -1.0 for _ in 1:nocc] ./ sqrt(nocc)
        v = Ψ * z
        exact = exact_action(Vxk, v)
        approx = ace_action(cached_op, v)
        ne = norm(exact)
        nd = norm(exact - approx)
        push!(vals, nd / max(ne, eps()))
        num += nd^2
        den += ne^2
    end
    estimate = sqrt(num / max(den, eps()))
    stderr = length(vals) > 1 ? std(vals) / sqrt(length(vals)) : 0.0
    score = estimate + 2 * stderr
    estimate, stderr, score
end


function approximate_exchange_energy(op, ψk, occk, maskk_occ)
    E = 0.0
    for i in maskk_occ
        v = @view ψk[:, i]
        E += 0.5 * occk[i] * real(dot(v, ace_action(op, v)))
    end
    E
end

function build_ace(basis, kpt, Vxk, ψk, occk, maskk_occ, ψk_occ_real, extra_count)
    nocc = length(maskk_occ)
    first_extra = nocc + 1
    last_extra = min(size(ψk, 2), nocc + max(0, extra_count))
    comp = first_extra <= last_extra ? collect(vcat(collect(maskk_occ), collect(first_extra:last_extra))) : collect(maskk_occ)
    ψcomp = @view ψk[:, comp]
    ψcomp_real = similar(ψk, basis.fft_size..., length(comp))
    for (ψn_real, idx) in zip(eachslice(ψcomp_real; dims=4), comp)
        DFTK.ifft!(ψn_real, basis, kpt, @view ψk[:, idx])
    end
    compressed = DFTK.compress_exchange(Vxk, ψcomp, ψcomp_real)
    E = 0.0
    for i in 1:nocc
        E += 0.5 * occk[i] * real(compressed.Mk[i, i])
    end
    compressed.opk, E, length(comp), Base.summarysize(compressed.opk.P) + Base.summarysize(compressed.opk.D)
end

mutable struct ErrorControlledAce <: DFTK.ExxAlgorithm
    defect_tol::Float64
    probe_count::Int
    max_reuse::Int
    min_extra::Int
    max_extra::Int
    extra_step::Int
    current_extra::Int
    call_count::Int
    last_rebuild::Int
    cached_op::Any
    calibrate::Bool
    calibration_stride::Int
    trace::Vector{Any}
end

function ErrorControlledAce(; defect_tol=3e-3, probe_count=1, max_reuse=5, min_extra=0, max_extra=8, extra_step=2, calibrate=false, calibration_stride=1)
    ErrorControlledAce(defect_tol, probe_count, max_reuse, min_extra, max_extra, extra_step, min_extra, 0, 0, nothing, calibrate, max(1, calibration_stride), Any[])
end

mutable struct PeriodicAce <: DFTK.ExxAlgorithm
    period::Int
    call_count::Int
    last_rebuild::Int
    cached_op::Any
    trace::Vector{Any}
end

PeriodicAce(period::Int) = PeriodicAce(max(1, period), 0, 0, nothing, Any[])

function DFTK.exx_operator(alg::PeriodicAce, basis::DFTK.PlaneWaveBasis{T}, kpt, interaction_kernel, ψk, occk, maskk_occ) where {T}
    t0 = time_ns()
    alg.call_count += 1
    Vxk, ψocc, ψocc_real, occocc = occupied_exchange_data(basis, kpt, interaction_kernel, ψk, occk, maskk_occ)
    rebuild = alg.cached_op === nothing || (alg.call_count - alg.last_rebuild >= alg.period)
    if rebuild
        opk, Ek, rank, opbytes = build_ace(basis, kpt, Vxk, ψk, occk, maskk_occ, ψocc_real, 0)
        alg.cached_op = opk
        alg.last_rebuild = alg.call_count
        action = "rebuild"
    else
        opk = alg.cached_op
        Ek = approximate_exchange_energy(opk, ψk, occk, maskk_occ)
        rank = size(opk.P, 2)
        opbytes = Base.summarysize(opk.P) + Base.summarysize(opk.D)
        action = "reuse"
    end
    elapsed = (time_ns() - t0) / 1e9
    push!(alg.trace, (call=alg.call_count, action=action, rank=rank, operator_bytes=opbytes, elapsed_s=elapsed))
    (; Ek, opk)
end

function DFTK.exx_operator(alg::ErrorControlledAce, basis::DFTK.PlaneWaveBasis{T}, kpt, interaction_kernel, ψk, occk, maskk_occ) where {T}
    t0 = time_ns()
    alg.call_count += 1
    Vxk, ψocc, ψocc_real, occocc = occupied_exchange_data(basis, kpt, interaction_kernel, ψk, occk, maskk_occ)
    nocc = length(maskk_occ)
    reuse_age = alg.cached_op === nothing ? 0 : alg.call_count - alg.last_rebuild
    forced_rebuild = alg.cached_op === nothing || reuse_age >= alg.max_reuse
    should_probe = alg.cached_op !== nothing && (!forced_rebuild || alg.calibrate)
    probe_seed = RNG_SEED + 104729 * alg.call_count
    if should_probe
        defect, defect_stderr, defect_score = stochastic_defect(Vxk, alg.cached_op, ψk, maskk_occ, alg.probe_count, probe_seed)
        probes_used = alg.probe_count
    else
        defect, defect_stderr, defect_score, probes_used = (NaN, NaN, NaN, 0)
    end
    defect_rebuild = alg.cached_op !== nothing && isfinite(defect_score) && defect_score > alg.defect_tol
    rebuild = forced_rebuild || defect_rebuild
    rebuild_reason = alg.cached_op === nothing ? "initial" : reuse_age >= alg.max_reuse ? "safety_cap" : defect_rebuild ? "defect" : "reuse"
    full_defect = NaN
    energy_defect = NaN
    replicate_defect = NaN
    calibration_probe3 = NaN
    calibration_probe3_stderr = NaN
    calibration_probe3_score = NaN
    calibrate_now = alg.calibrate && alg.cached_op !== nothing && ((alg.call_count - 1) % alg.calibration_stride == 0)
    if calibrate_now
        full_defect = full_occupied_defect(Vxk, alg.cached_op, ψk, maskk_occ)
        calibration_probe3, calibration_probe3_stderr, calibration_probe3_score = stochastic_defect(Vxk, alg.cached_op, ψk, maskk_occ, 3, probe_seed)
        replicate_defect = stochastic_defect(Vxk, alg.cached_op, ψk, maskk_occ, 1, RNG_SEED + 999983 + 104729 * alg.call_count)[1]
        Eexact = DFTK.exx_energy_only(basis, kpt, interaction_kernel, ψocc_real, occocc)
        Eapprox = approximate_exchange_energy(alg.cached_op, ψk, occk, maskk_occ)
        energy_defect = abs(Eapprox - Eexact)
    end
    rank = nocc + alg.current_extra
    opbytes = 0
    if rebuild
        if isfinite(defect_score)
            if defect_score > 3 * alg.defect_tol
                alg.current_extra = min(alg.max_extra, alg.current_extra + alg.extra_step)
            elseif defect_score < alg.defect_tol / 4
                alg.current_extra = max(alg.min_extra, alg.current_extra - alg.extra_step)
            end
        end
        opk, Ek, rank, opbytes = build_ace(basis, kpt, Vxk, ψk, occk, maskk_occ, ψocc_real, alg.current_extra)
        alg.cached_op = opk
        alg.last_rebuild = alg.call_count
        action = "rebuild"
    else
        opk = alg.cached_op
        Ek = approximate_exchange_energy(opk, ψk, occk, maskk_occ)
        rank = size(opk.P, 2)
        opbytes = Base.summarysize(opk.P) + Base.summarysize(opk.D)
        action = "reuse"
    end
    elapsed = (time_ns() - t0) / 1e9
    push!(alg.trace, (call=alg.call_count, action=action, probe_defect=defect, probe_stderr=defect_stderr, defect_score=defect_score, full_defect=full_defect, replicate_defect=replicate_defect, calibration_probe3=calibration_probe3, calibration_probe3_stderr=calibration_probe3_stderr, calibration_probe3_score=calibration_probe3_score, energy_defect=energy_defect, rank=rank, operator_bytes=opbytes, elapsed_s=elapsed, reuse_age=reuse_age, rebuild_reason=rebuild_reason, probes_used=probes_used))
    (; Ek, opk)
end

function model_for(cfg; hybrid=true)
    lattice, atoms, positions = cfg.geom
    fun = cfg.functional == "HSE06" ? HSE() : PBE0()
    if hybrid
        model_DFT(lattice, atoms, positions; functionals=fun)
    else
        model_DFT(lattice, atoms, positions; functionals=PBE())
    end
end

function prepare_pbe(cfg)
    key = (cfg.name, cfg.Ecut)
    haskey(PBE_CACHE, key) && return PBE_CACHE[key]
    model = model_for(cfg; hybrid=false)
    basis = PlaneWaveBasis(model; Ecut=cfg.Ecut, kgrid=[1,1,1], architecture=DFTK.CPU())
    rows = Any[]
    cb = info -> begin
        info.stage == :iterate || return
        push!(rows, (length(rows)+1, isempty(info.history_Δρ) ? NaN : info.history_Δρ[end], info.runtime_ns/1e9))
    end
    t0 = time_ns()
    scfres = self_consistent_field(basis; tol=1e-8, maxiter=80, callback=cb, seed=DFTK_SEED)
    wall = (time_ns()-t0)/1e9
    PBE_CACHE[key] = (scfres, wall)
    PBE_CACHE[key]
end

function run_hybrid(cfg, pberes, method, exxalg)
    model = model_for(cfg; hybrid=true)
    basis = PlaneWaveBasis(model; Ecut=cfg.Ecut, kgrid=[1,1,1], architecture=DFTK.CPU())
    hist = Any[]
    cb = info -> begin
        info.stage == :iterate || return
        push!(hist, (cfg.name, method, length(hist)+1, info.n_iter, isempty(info.history_Δρ) ? NaN : info.history_Δρ[end], info.runtime_ns/1e9))
    end
    t0 = time_ns()
    scfres = self_consistent_field(basis;
        tol=1e-6,
        maxiter=80,
        damping=0.30,
        solver=ScfDampingSolver(),
        mixing=LdosMixing(),
        callback=cb,
        seed=DFTK_SEED,
        compute_consistent_energies=false,
        exxalg=exxalg,
        ρ=pberes.ρ,
        ψ=pberes.ψ,
        occupation=pberes.occupation,
        eigenvalues=pberes.eigenvalues,
    )
    wall = (time_ns()-t0)/1e9
    inner = isempty(hist) ? NaN : hist[end][6]
    final_residual = isempty(hist) ? NaN : hist[end][5]
    summary = (cfg.name, cfg.functional, method, scfres.converged, length(hist), scfres.n_iter, scfres.n_matvec, inner, wall, scfres.energies.total, final_residual)
    hist, summary
end

function run_main_suite(cfgs)
    hist_rows = Any[]
    summary_rows = Any[]
    trace_rows = Any[]
    system_rows = Any[]
    for cfg in cfgs
        append_log("system_start", cfg.name, cfg.functional, nowstamp())
        pberes, pbewall = prepare_pbe(cfg)
        lattice, atoms, positions = cfg.geom
        push!(system_rows, (cfg.name, cfg.functional, length(atoms), cfg.Ecut, "1x1x1", pbewall))
        methods = [
            ("Vanilla_EXX", VanillaExx()),
            ("ACE_occupied", AceExx(sketch_with_extra_orbitals=false)),
            ("ACE_full", AceExx(sketch_with_extra_orbitals=true)),
            ("AdaptiveEXX", ErrorControlledAce(defect_tol=3e-3, probe_count=1, max_reuse=5, min_extra=0, max_extra=8)),
        ]
        for (name, alg) in methods
            GC.gc()
            hist, summary = run_hybrid(cfg, pberes, name, alg)
            append!(hist_rows, hist)
            push!(summary_rows, summary)
            if alg isa ErrorControlledAce
                for r in alg.trace
                    push!(trace_rows, (cfg.name, cfg.functional, r.call, r.action, r.probe_defect, r.probe_stderr, r.defect_score, r.full_defect, r.replicate_defect, r.energy_defect, r.rank, r.operator_bytes/2.0^20, r.elapsed_s, r.probes_used))
                end
            end
            append_log("method_done", cfg.name, name, "converged", summary[4], "evaluations", summary[5], "inner_s", summary[8], "wall_s", summary[9], "energy_Ha", summary[10])
        end
    end
    hist_rows, summary_rows, trace_rows, system_rows
end

function run_calibration(cfgs)
    rows = Any[]
    for cfg in cfgs
        pberes, _ = prepare_pbe(cfg)
        stride = cfg.name == "Si16" ? 2 : 1
        alg = ErrorControlledAce(defect_tol=3e-3, probe_count=1, max_reuse=5, min_extra=0, max_extra=6, calibrate=true, calibration_stride=stride)
        run_hybrid(cfg, pberes, "AdaptiveEXX_calibration", alg)
        for r in alg.trace
            isfinite(r.full_defect) || continue
            probe1_rel = r.full_defect > 0 ? abs(r.probe_defect - r.full_defect) / r.full_defect : NaN
            probe3_rel = isfinite(r.calibration_probe3) && r.full_defect > 0 ? abs(r.calibration_probe3 - r.full_defect) / r.full_defect : NaN
            replicate_rel = isfinite(r.replicate_defect) && r.full_defect > 0 ? abs(r.replicate_defect - r.full_defect) / r.full_defect : NaN
            push!(rows, (cfg.name, r.call, r.action, r.probe_defect, probe1_rel, r.calibration_probe3, r.calibration_probe3_stderr, r.calibration_probe3_score, probe3_rel, r.full_defect, r.replicate_defect, replicate_rel, r.energy_defect*HARTREE2EV*1000, r.rank, r.elapsed_s))
        end
    end
    rows
end

function run_sensitivity(cfg, reference_energy; settings=nothing)
    pberes, _ = prepare_pbe(cfg)
    rows = Any[]
    settings === nothing && (settings = [(1e-2,3,5), (3e-3,3,5), (1e-3,3,5), (3e-3,1,5), (3e-3,5,5)])
    for (tol, probes, reuse) in settings
        alg = ErrorControlledAce(defect_tol=tol, probe_count=probes, max_reuse=reuse, min_extra=0, max_extra=8)
        hist, summary = run_hybrid(cfg, pberes, "adaptive_tol$(tol)_probe$(probes)", alg)
        nreb = count(r -> r.action == "rebuild", alg.trace)
        nreuse = count(r -> r.action == "reuse", alg.trace)
        maxdef = isempty(alg.trace) ? NaN : maximum(r -> isfinite(r.defect_score) ? r.defect_score : 0.0, alg.trace)
        used_probes = [Float64(r.probes_used) for r in alg.trace if r.probes_used > 0]
        mean_used_probes = isempty(used_probes) ? 0.0 : mean(used_probes)
        total_used_probes = sum(used_probes)
        nat = length(cfg.geom[2])
        eerr = abs(summary[10] - reference_energy) * HARTREE2EV * 1000 / nat
        push!(rows, (cfg.name, tol, probes, reuse, summary[4], summary[5], summary[8], summary[9], summary[10], eerr, nreb, nreuse, maxdef, mean_used_probes, total_used_probes))
        append_log("probe_sensitivity_done", cfg.name, "probes", probes, "tol", tol, "converged", summary[4], "evaluations", summary[5], "wall_s", summary[9], "monitoring_probes", total_used_probes)
    end
    rows
end

function run_periodic_baselines(cfg, reference_energy)
    pberes, _ = prepare_pbe(cfg)
    rows = Any[]
    nat = length(cfg.geom[2])
    for period in (2, 3, 5)
        alg = PeriodicAce(period)
        _, summary = run_hybrid(cfg, pberes, "PeriodicACE_$(period)", alg)
        rebuilds = count(r -> r.action == "rebuild", alg.trace)
        reuses = count(r -> r.action == "reuse", alg.trace)
        eerr = abs(summary[10] - reference_energy) * HARTREE2EV * 1000 / nat
        push!(rows, (cfg.name, "PeriodicACE_$(period)", period, summary[4], summary[5], summary[8], summary[9], summary[10], eerr, rebuilds, reuses, sum(r.elapsed_s for r in alg.trace)))
    end
    rows
end

function run_size_scaling(main_summary, cfg16)
    rows = Any[]
    for system in ("Si2", "Si8")
        for method in ("ACE_occupied", "AdaptiveEXX")
            hit = first(r for r in main_summary if r[1] == system && r[3] == method)
            atoms = system == "Si2" ? 2 : 8
            push!(rows, (system, atoms, method, hit[4], hit[5], hit[8], hit[9], hit[10]))
        end
    end
    pberes, _ = prepare_pbe(cfg16)
    for (name, alg) in (("ACE_occupied", AceExx(sketch_with_extra_orbitals=false)), ("AdaptiveEXX", ErrorControlledAce(defect_tol=3e-3, probe_count=1, max_reuse=5, min_extra=0, max_extra=8)))
        _, summary = run_hybrid(cfg16, pberes, name, alg)
        push!(rows, (cfg16.name, length(cfg16.geom[2]), name, summary[4], summary[5], summary[8], summary[9], summary[10]))
        append_log("scaling_done", cfg16.name, name, "converged", summary[4], "wall_s", summary[9])
    end
    rows
end

function warmup()
    cfg = (name="warmup", geom=diamond(:Si, 5.431), Ecut=12.0, functional="PBE0")
    pmodel = model_for(cfg; hybrid=false)
    pbasis = PlaneWaveBasis(pmodel; Ecut=cfg.Ecut, kgrid=[1,1,1], architecture=DFTK.CPU())
    pberes = self_consistent_field(pbasis; tol=1e-2, maxiter=2, callback=info -> nothing, seed=DFTK_SEED)
    hmodel = model_for(cfg; hybrid=true)
    hbasis = PlaneWaveBasis(hmodel; Ecut=cfg.Ecut, kgrid=[1,1,1], architecture=DFTK.CPU())
    for alg in (VanillaExx(), AceExx(sketch_with_extra_orbitals=false), ErrorControlledAce(defect_tol=1e-2, probe_count=1, max_reuse=2))
        try
            self_consistent_field(hbasis; tol=1e-2, maxiter=2, damping=0.3, solver=ScfDampingSolver(), mixing=LdosMixing(), callback=info -> nothing, seed=DFTK_SEED, compute_consistent_energies=false, exxalg=alg, ρ=pberes.ρ, ψ=pberes.ψ, occupation=pberes.occupation, eigenvalues=pberes.eigenvalues)
        catch
        end
    end
    GC.gc()
end

function main()
    println("Precompiling DFTK and exact-exchange kernels; warmup is excluded from reported timings.")
    warmup()
    write_log_header()
    tmain = time_ns()
    append_log("event", "start", nowstamp())
    main_hist, main_summary, main_trace, system_rows = run_main_suite(systems())
    hse_hist, hse_summary, hse_trace, hse_system_rows = run_main_suite(hse_systems())
    reference = Dict{Tuple{String,String},Float64}()
    for r in main_summary
        r[3] == "Vanilla_EXX" && (reference[(r[1],r[2])] = r[10])
    end
    error_rows = Any[]
    for r in main_summary
        e0 = reference[(r[1],r[2])]
        nat = first(x[3] for x in system_rows if x[1] == r[1])
        push!(error_rows, (r[1], r[3], abs(r[10]-e0)*HARTREE2EV*1000/nat, r[8], r[9], r[4]))
    end
    hse_reference = Dict{Tuple{String,String},Float64}()
    for r in hse_summary
        r[3] == "Vanilla_EXX" && (hse_reference[(r[1],r[2])] = r[10])
    end
    hse_error_rows = Any[]
    for r in hse_summary
        e0 = hse_reference[(r[1],r[2])]
        nat = first(x[3] for x in hse_system_rows if x[1] == r[1])
        push!(hse_error_rows, (r[1], r[3], abs(r[10]-e0)*HARTREE2EV*1000/nat, r[8], r[9], r[4]))
    end
    cfg16 = si16_system()
    calibration_rows = run_calibration([systems()[1], systems()[3], cfg16])
    si8_ref = reference[("Si8", "PBE0")]
    mgo2_ref = reference[("MgO2", "PBE0")]
    sensitivity_rows = run_sensitivity(systems()[5], si8_ref)
    append!(sensitivity_rows, run_sensitivity(systems()[3], mgo2_ref; settings=[(3e-3,1,5), (3e-3,3,5), (3e-3,5,5)]))
    periodic_rows = Any[]
    for cfg in (systems()[3], systems()[5])
        refE = reference[(cfg.name, "PBE0")]
        append!(periodic_rows, run_periodic_baselines(cfg, refE))
        nat = length(cfg.geom[2])
        for method in ("ACE_occupied", "AdaptiveEXX")
            hit = first(r for r in main_summary if r[1] == cfg.name && r[3] == method)
            eerr = abs(hit[10] - refE) * HARTREE2EV * 1000 / nat
            push!(periodic_rows, (cfg.name, method, 0, hit[4], hit[5], hit[8], hit[9], hit[10], eerr, NaN, NaN, NaN))
        end
    end
    scaling_rows = run_size_scaling(main_summary, cfg16)
    controller_rows = Any[]
    for system in sort(unique([r[1] for r in main_trace]))
        rr = [r for r in main_trace if r[1] == system]
        isempty(rr) && continue
        rebuilds = count(r -> r[4] == "rebuild", rr)
        reuses = count(r -> r[4] == "reuse", rr)
        ranks = [Float64(r[11]) for r in rr if isfinite(Float64(r[11]))]
        memories = [Float64(r[12]) for r in rr if isfinite(Float64(r[12]))]
        times = [Float64(r[13]) for r in rr if isfinite(Float64(r[13]))]
        scores = [Float64(r[7]) for r in rr if isfinite(Float64(r[7]))]
        probes = [Float64(r[14]) for r in rr if isfinite(Float64(r[14])) && Float64(r[14]) > 0]
        push!(controller_rows, (system, rr[1][2], length(rr), rebuilds, reuses, isempty(ranks) ? NaN : mean(ranks), isempty(memories) ? NaN : maximum(memories), sum(times), isempty(scores) ? NaN : maximum(scores), isempty(probes) ? 0.0 : mean(probes), sum(probes)))
    end
    probe_efficiency_rows = Any[]
    for system in ("Si8", "MgO2")
        for probes in (1, 3, 5)
            row = first(r for r in sensitivity_rows if r[1] == system && abs(r[2]-3e-3) < 1e-12 && r[3] == probes && r[4] == 5)
            push!(probe_efficiency_rows, (system, "Fixed $(probes)", row[14], row[15], row[6], row[8], row[10], row[11], row[12]))
        end
    end
    write_tsv(joinpath(MAIN_DIR, "fig1_defect_calibration.dat"),
        ["system","call","action","probe1_defect","probe1_relative_error","probe3_defect","probe3_stderr","probe3_score","probe3_relative_error","full_occupied_defect","independent_probe1_defect","independent_probe1_relative_error","exchange_energy_defect_meV","ace_rank","exx_call_time_s"], calibration_rows)
    write_tsv(joinpath(MAIN_DIR, "fig2_scf_convergence.dat"),
        ["system","method","evaluation","accepted_iteration","density_residual","elapsed_s"], main_hist)
    write_tsv(joinpath(MAIN_DIR, "fig3_cost_accuracy.dat"),
        ["system","method","energy_deviation_from_vanilla_meV_per_atom","inner_scf_time_s","total_wall_time_s","converged"], error_rows)
    write_tsv(joinpath(MAIN_DIR, "fig4_size_scaling.dat"),
        ["system","atoms","method","converged","scf_evaluations","inner_scf_time_s","total_wall_time_s","final_energy_Ha"], scaling_rows)
    write_tsv(joinpath(MAIN_DIR, "fig4_probe_efficiency.dat"),
        ["system","strategy","mean_probes_when_monitored","total_monitoring_probes","scf_evaluations","total_wall_time_s","energy_deviation_from_vanilla_meV_per_atom","rebuilds","reuses"], probe_efficiency_rows)
    write_tsv(joinpath(MAIN_DIR, "fig5_controller_efficiency.dat"),
        ["system","functional","calls","rebuilds","reuses","mean_ace_rank","max_operator_memory_MiB","total_exx_call_time_s","max_defect_score","mean_probes_when_monitored","total_monitoring_probes"], controller_rows)
    write_tsv(joinpath(MAIN_DIR, "fig6_hse_transfer.dat"),
        ["system","method","energy_deviation_from_vanilla_meV_per_atom","inner_scf_time_s","total_wall_time_s","converged"], hse_error_rows)
    write_tsv(joinpath(MAIN_DIR, "table1_systems.dat"),
        ["system","functional","atoms","Ecut_Ha","kgrid","pbe_initialization_time_s"], system_rows)
    write_tsv(joinpath(MAIN_DIR, "table2_performance.dat"),
        ["system","functional","method","converged","scf_evaluations","accepted_iterations","hamiltonian_matvecs","inner_scf_time_s","total_wall_time_s","final_energy_Ha","final_density_residual"], main_summary)
    write_tsv(joinpath(MAIN_DIR, "table3_controller_statistics.dat"),
        ["system","functional","calls","rebuilds","reuses","mean_ace_rank","max_operator_memory_MiB","total_exx_call_time_s","max_defect_score","mean_probes_when_monitored","total_monitoring_probes"], controller_rows)
    write_tsv(joinpath(SI_DIR, "figS1_hse_convergence.dat"),
        ["system","method","evaluation","accepted_iteration","density_residual","elapsed_s"], hse_hist)
    write_tsv(joinpath(SI_DIR, "tableS1_hse_performance.dat"),
        ["system","functional","method","converged","scf_evaluations","accepted_iterations","hamiltonian_matvecs","inner_scf_time_s","total_wall_time_s","final_energy_Ha","final_density_residual"], hse_summary)
    write_tsv(joinpath(SI_DIR, "figS2_sensitivity.dat"),
        ["system","defect_tolerance","probe_count","max_reuse","converged","scf_evaluations","inner_scf_time_s","total_wall_time_s","final_energy_Ha","energy_deviation_from_vanilla_meV_per_atom","rebuilds","reuses","max_defect_score","mean_probes_when_monitored","total_monitoring_probes"], sensitivity_rows)
    write_tsv(joinpath(SI_DIR, "tableS2_sensitivity.dat"),
        ["system","defect_tolerance","probe_count","max_reuse","converged","scf_evaluations","inner_scf_time_s","total_wall_time_s","final_energy_Ha","energy_deviation_from_vanilla_meV_per_atom","rebuilds","reuses","max_defect_score","mean_probes_when_monitored","total_monitoring_probes"], sensitivity_rows)
    write_tsv(joinpath(SI_DIR, "figS3_hse_adaptive_trace.dat"),
        ["system","functional","call","action","probe_defect","probe_stderr","defect_score","full_defect","replicate_defect","energy_defect_Ha","ace_rank","operator_memory_MiB","exx_call_time_s","probes_used"], hse_trace)
    write_tsv(joinpath(SI_DIR, "figS4_all_pbe_convergence.dat"),
        ["system","method","evaluation","accepted_iteration","density_residual","elapsed_s"], main_hist)
    write_tsv(joinpath(SI_DIR, "figS5_pbe_energy_error.dat"),
        ["system","method","energy_deviation_from_vanilla_meV_per_atom","inner_scf_time_s","total_wall_time_s","converged"], error_rows)
    write_tsv(joinpath(SI_DIR, "figS6_method_speedup.dat"),
        ["system","method","energy_deviation_from_vanilla_meV_per_atom","inner_scf_time_s","total_wall_time_s","converged"], error_rows)
    write_tsv(joinpath(SI_DIR, "figS7_calibration_relative_error.dat"),
        ["system","call","action","probe1_defect","probe1_relative_error","probe3_defect","probe3_stderr","probe3_score","probe3_relative_error","full_occupied_defect","independent_probe1_defect","independent_probe1_relative_error","exchange_energy_defect_meV","ace_rank","exx_call_time_s"], calibration_rows)
    write_tsv(joinpath(SI_DIR, "figS8_sensitivity_actions.dat"),
        ["system","defect_tolerance","probe_count","max_reuse","converged","scf_evaluations","inner_scf_time_s","total_wall_time_s","final_energy_Ha","energy_deviation_from_vanilla_meV_per_atom","rebuilds","reuses","max_defect_score","mean_probes_when_monitored","total_monitoring_probes"], sensitivity_rows)
    write_tsv(joinpath(SI_DIR, "tableS3_hse_accuracy.dat"),
        ["system","method","energy_deviation_from_vanilla_meV_per_atom","inner_scf_time_s","total_wall_time_s","converged"], hse_error_rows)
    write_tsv(joinpath(SI_DIR, "tableS4_estimator_calibration.dat"),
        ["system","call","action","probe1_defect","probe1_relative_error","probe3_defect","probe3_stderr","probe3_score","probe3_relative_error","full_occupied_defect","independent_probe1_defect","independent_probe1_relative_error","exchange_energy_defect_meV","ace_rank","exx_call_time_s"], calibration_rows)
    write_tsv(joinpath(SI_DIR, "figS9_periodic_rebuild.dat"),
        ["system","method","period","converged","scf_evaluations","inner_scf_time_s","total_wall_time_s","final_energy_Ha","energy_deviation_from_vanilla_meV_per_atom","rebuilds","reuses","total_exx_call_time_s"], periodic_rows)
    write_tsv(joinpath(SI_DIR, "tableS5_periodic_rebuild.dat"),
        ["system","method","period","converged","scf_evaluations","inner_scf_time_s","total_wall_time_s","final_energy_Ha","energy_deviation_from_vanilla_meV_per_atom","rebuilds","reuses","total_exx_call_time_s"], periodic_rows)
    write_tsv(joinpath(SI_DIR, "figS10_adaptive_trace.dat"),
        ["system","functional","call","action","probe_defect","probe_stderr","defect_score","full_defect","replicate_defect","energy_defect_Ha","ace_rank","operator_memory_MiB","exx_call_time_s","probes_used"], main_trace)
    append_log("event", "finish", nowstamp(), "production_wall_s", (time_ns()-tmain)/1e9)
    println("AdaptiveEXX.jl completed.")
    println("Run: python redraw_paper.py")
end

main()
