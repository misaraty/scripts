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
const AUTO_INSTALL = lowercase(get(ENV, "GOALDFT_AUTO_INSTALL", "1")) in ("1", "true", "yes", "on")
const PSP_NAME = "dojo.nc.sr.pbe.v0_4_1.standard.upf"
const ANG2BOHR = 1.8897261254578281
const HARTREE2EV = 27.211386245988
const BOHR2ANG = 0.529177210903
const HA2MEV = HARTREE2EV * 1000
const HAPERBOHR2EVPERANG = HARTREE2EV / BOHR2ANG
const MAIN_DIR = joinpath(ROOT, "data", "main")
const SI_DIR = joinpath(ROOT, "data", "si")
const FIG_DIR = joinpath(ROOT, "figures")
const LOG_FILE = joinpath(ROOT, "log.dat")
const RNG_SEED = 20260811
const DFTK_SEED = UInt64(RNG_SEED)
const ENERGY_TARGETS_MEV = [10.0, 5.0, 2.0]
const FORCE_TARGETS_EVANG = [0.05, 0.02]
const SAFETY_FACTOR = 0.65
const NEAR_THRESHOLD_FACTOR = 1.10
const NEAR_THRESHOLD_VERIFY_INFLATION = 1.10
const MAX_CONTROL_STEPS = 10
const OCCUPATION_WARMSTART_TEMPERATURE_HA = 1e-3
const OCCUPATION_WARMSTART_MAXITER = 100

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

function nowstamp()
    Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS")
end

function write_log_header()
    open(LOG_FILE, "w") do io
        println(io, "project\tGoalDFT.jl")
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
        println(io, "total_memory_GiB\t$(round(Sys.total_memory()/2.0^30, digits=3))")
        println(io, "free_memory_GiB\t$(round(Sys.free_memory()/2.0^30, digits=3))")
        println(io, "rng_seed\t$RNG_SEED")
        println(io, "dftk_seed_type\t$(typeof(DFTK_SEED))")
        println(io, "safety_factor\t$SAFETY_FACTOR")
        println(io, "near_threshold_factor\t$NEAR_THRESHOLD_FACTOR")
        println(io, "near_threshold_verify_inflation\t$NEAR_THRESHOLD_VERIFY_INFLATION")
        println(io, "controller_policy\tlazy k-point probing with decay prediction, near-threshold direct verification, boundary stop verification, one anchored SCF probe, shared batch cache")
        println(io, "production_temperature_Ha\t0.0")
        println(io, "occupation_warmstart_temperature_Ha\t$OCCUPATION_WARMSTART_TEMPERATURE_HA")
        println(io, "occupation_warmstart_policy\tfinite-temperature density preconditioner only after zero-temperature occupation failure; final reported SCF remains zero-temperature")
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

function diamond(sym, a_ang; displacement=nothing)
    a = a_ang * ANG2BOHR
    lattice = (a / 2) .* [0.0 1.0 1.0; 1.0 0.0 1.0; 1.0 1.0 0.0]
    X = ElementPsp(sym, PSP)
    positions = [ones(3)/8, -ones(3)/8]
    if displacement !== nothing
        positions[1] = positions[1] .+ displacement
    end
    lattice, [X, X], positions
end

function rocksalt(a_sym, b_sym, a_ang; displacement=nothing)
    a = a_ang * ANG2BOHR
    lattice = (a / 2) .* [0.0 1.0 1.0; 1.0 0.0 1.0; 1.0 1.0 0.0]
    A = ElementPsp(a_sym, PSP)
    B = ElementPsp(b_sym, PSP)
    positions = [zeros(3), fill(0.5, 3)]
    if displacement !== nothing
        positions[2] = positions[2] .+ displacement
    end
    lattice, [A, B], positions
end

function tio2()
    Ti = ElementPsp(:Ti, PSP)
    O = ElementPsp(:O, PSP)
    atoms = [Ti, Ti, O, O, O, O]
    positions = [
        [0.5, 0.5, 0.5],
        [0.0, 0.0, 0.0],
        [0.19542, 0.80458, 0.5],
        [0.80458, 0.19542, 0.5],
        [0.30458, 0.30458, 0.0],
        [0.69542, 0.69542, 0.0],
    ]
    positions[1] = positions[1] .+ [-0.022, 0.028, 0.035]
    lattice = [8.79341 0.0 0.0; 0.0 8.79341 0.0; 0.0 0.0 5.61098]
    lattice, atoms, positions
end

function energy_systems()
    [
        (name="Si", geom=diamond(:Si, 5.431), start_Ecut=28.0, ref_Ecut=55.0, truth_Ecut=70.0, start_nk=2, ref_nk=8, truth_nk=10, start_tol=1e-4),
        (name="C", geom=diamond(:C, 3.567), start_Ecut=32.0, ref_Ecut=65.0, truth_Ecut=80.0, start_nk=2, ref_nk=8, truth_nk=10, start_tol=1e-4),
        (name="MgO", geom=rocksalt(:Mg, :O, 4.212), start_Ecut=32.0, ref_Ecut=65.0, truth_Ecut=80.0, start_nk=2, ref_nk=8, truth_nk=10, start_tol=1e-4),
        (name="LiF", geom=rocksalt(:Li, :F, 4.027), start_Ecut=36.0, ref_Ecut=70.0, truth_Ecut=85.0, start_nk=2, ref_nk=8, truth_nk=10, start_tol=1e-4),
    ]
end

function force_systems()
    [
        (name="Si_displaced", geom=diamond(:Si, 5.431; displacement=[0.012,-0.009,0.007]), start_Ecut=30.0, ref_Ecut=60.0, truth_Ecut=75.0, start_nk=2, ref_nk=8, truth_nk=10, start_tol=1e-5),
        (name="MgO_displaced", geom=rocksalt(:Mg, :O, 4.212; displacement=[0.015,-0.010,0.008]), start_Ecut=34.0, ref_Ecut=68.0, truth_Ecut=82.0, start_nk=2, ref_nk=8, truth_nk=10, start_tol=1e-5),
        (name="TiO2_displaced", geom=tio2(), start_Ecut=32.0, ref_Ecut=65.0, truth_Ecut=80.0, start_nk=2, ref_nk=6, truth_nk=8, start_tol=1e-5),
    ]
end

function build_model(cfg; temperature=0.0)
    lattice, atoms, positions = cfg.geom
    if temperature > 0
        model_DFT(lattice, atoms, positions; functionals=PBE(), temperature=temperature, smearing=DFTK.Smearing.FermiDirac())
    else
        model_DFT(lattice, atoms, positions; functionals=PBE())
    end
end

function is_occupation_error(err)
    msg = sprint(showerror, err)
    occursin("Unable to find non-fractional occupations", msg) || occursin("cannot be attained by filling states", msg)
end

function run_scf(cfg, Ecut, nk, tol, goal)
    model = build_model(cfg)
    basis = PlaneWaveBasis(model; Ecut=Ecut, kgrid=[nk,nk,nk], architecture=DFTK.CPU())
    histE = Float64[]
    histR = Float64[]
    cb = info -> begin
        info.stage == :iterate || return
        hasproperty(info, :energies) && push!(histE, info.energies.total)
        push!(histR, isempty(info.history_Δρ) ? NaN : info.history_Δρ[end])
    end
    GC.gc()
    t0 = time_ns()
    used_warmstart = false
    warmstart_wall = 0.0
    warmstart_converged = true
    scfres = try
        self_consistent_field(basis; tol=tol, maxiter=120, callback=cb, seed=DFTK_SEED)
    catch err
        is_occupation_error(err) || rethrow()
        used_warmstart = true
        empty!(histE)
        empty!(histR)
        append_log("occupation_warmstart_start", cfg.name, goal, "Ecut", Ecut, "kgrid", nk, "target_tol", tol, "temperature_Ha", OCCUPATION_WARMSTART_TEMPERATURE_HA)
        warm_model = build_model(cfg; temperature=OCCUPATION_WARMSTART_TEMPERATURE_HA)
        warm_basis = PlaneWaveBasis(warm_model; Ecut=Ecut, kgrid=[nk,nk,nk], architecture=DFTK.CPU())
        tw = time_ns()
        warm_tol = max(tol, 1e-6)
        warm = self_consistent_field(warm_basis; tol=warm_tol, maxiter=OCCUPATION_WARMSTART_MAXITER, callback=identity, seed=DFTK_SEED)
        warmstart_wall = (time_ns()-tw)/1e9
        warmstart_converged = warm.converged
        append_log("occupation_warmstart_done", cfg.name, goal, "converged", warm.converged, "wall_s", warmstart_wall, "warm_tol", warm_tol)
        self_consistent_field(basis; ρ=warm.ρ, tol=tol, maxiter=120, callback=cb, seed=DFTK_SEED)
    end
    forces = goal == :force ? compute_forces_cart(scfres) : nothing
    wall = (time_ns()-t0)/1e9
    append_log("scf_run", cfg.name, goal, "Ecut", Ecut, "kgrid", nk, "tol", tol, "converged", scfres.converged, "wall_s", wall, "occupation_warmstart", used_warmstart, "warmstart_wall_s", warmstart_wall, "warmstart_converged", warmstart_converged)
    (; scfres, wall, Ecut, nk, tol, goal, histE, histR, energy=scfres.energies.total, forces, used_warmstart, warmstart_wall, warmstart_converged)
end

function max_force_difference(a, b)
    maximum(norm(a[i] - b[i]) for i in eachindex(a))
end

function property_difference(a, b, cfg, goal)
    if goal == :energy
        abs(a.energy - b.energy) * HA2MEV / length(cfg.geom[2])
    else
        max_force_difference(a.forces, b.forces) * HAPERBOHR2EVPERANG
    end
end

function run_key(Ecut, nk, tol, goal)
    (round(Ecut, digits=8), nk, round(tol, sigdigits=4), goal)
end

function get_run!(cache, used, spent, cfg, Ecut, nk, tol, goal)
    key = run_key(Ecut, nk, tol, goal)
    if !haskey(cache, key)
        cache[key] = run_scf(cfg, Ecut, nk, tol, goal)
    end
    if !(key in used)
        push!(used, key)
        spent[] += cache[key].wall
    end
    cache[key]
end

function convert_refined_force_delta(model, dF)
    DFTK.covector_red_to_cart.(model, dF)
end

function basis_estimate!(refcache, used_ref, spent, cfg, current, new_Ecut, goal)
    key = (run_key(current.Ecut, current.nk, current.tol, goal), round(new_Ecut, digits=8))
    if !haskey(refcache, key)
        model = current.scfres.basis.model
        basis_ref = PlaneWaveBasis(model; Ecut=new_Ecut, kgrid=[current.nk,current.nk,current.nk], architecture=DFTK.CPU())
        t0 = time_ns()
        refinement = refine_scfres(current.scfres, basis_ref; tol=min(1e-7, current.tol/10), callback=identity)
        if goal == :energy
            d = refine_energies(refinement).dE.total
            estimate = abs(d) * HA2MEV / length(cfg.geom[2])
        else
            dF = refine_forces(refinement).dF
            dF_cart = convert_refined_force_delta(basis_ref.model, dF)
            estimate = maximum(norm.(dF_cart)) * HAPERBOHR2EVPERANG
        end
        refcache[key] = (; estimate, wall=(time_ns()-t0)/1e9)
    end
    if !(key in used_ref)
        push!(used_ref, key)
        spent[] += refcache[key].wall
    end
    refcache[key]
end

function reference_run(cfg, goal)
    run_scf(cfg, cfg.ref_Ecut, cfg.ref_nk, 1e-9, goal)
end

function truth_run(cfg, goal)
    run_scf(cfg, cfg.truth_Ecut, cfg.truth_nk, 1e-10, goal)
end

function next_ecut(cfg, Ecut)
    min(cfg.ref_Ecut, max(Ecut + 5.0, 1.18 * Ecut))
end

function next_tol(tol)
    max(1e-9, tol / 10)
end

function choose_action(basis_est, k_est, scf_est, basis_cost, k_cost, scf_cost, Ecut, new_Ecut, nk, max_nk, tol)
    items = [
        (:Ecut, basis_est / max(basis_cost, 1e-9), basis_est, new_Ecut > Ecut + 1e-10),
        (:kgrid, k_est / max(k_cost, 1e-9), k_est, nk < max_nk),
        (:scf, scf_est / max(scf_cost, 1e-9), scf_est, tol > 1.01e-9),
    ]
    valid = filter(x -> x[4] && isfinite(x[2]), items)
    isempty(valid) && return :limit
    sort!(valid; by=x -> (x[2], x[3]), rev=true)
    valid[1][1]
end

function predict_next_k_error(last_change, old_nk, new_nk)
    isfinite(last_change) || return Inf
    last_change * (old_nk / new_nk)^3
end

function adaptive_goal(cfg, target, goal, truth, stringent_wall, run_cache, refcache)
    used = Set{Any}()
    used_ref = Set{Any}()
    spent = Ref(0.0)
    Ecut = cfg.start_Ecut
    nk = cfg.start_nk
    tol = cfg.start_tol
    max_nk = cfg.ref_nk + 1
    current = get_run!(run_cache, used, spent, cfg, Ecut, nk, tol, goal)
    tol_anchor = next_tol(tol)
    tight_anchor = tol_anchor < tol ? get_run!(run_cache, used, spent, cfg, Ecut, nk, tol_anchor, goal) : current
    scf_anchor = tol_anchor < tol ? property_difference(current, tight_anchor, cfg, goal) : 0.0
    k_est = Inf
    k_source = "uninitialized"
    k_candidate = nothing
    k_candidate_n = nk
    if nk < cfg.ref_nk
        k_candidate_n = nk + 1
        k_candidate = get_run!(run_cache, used, spent, cfg, Ecut, k_candidate_n, tol, goal)
        k_est = property_difference(current, k_candidate, cfg, goal)
        k_source = "direct_probe"
    else
        k_est = 0.0
        k_source = "at_limit"
    end
    trace = Any[]
    for step in 1:MAX_CONTROL_STEPS
        actual = property_difference(current, truth, cfg, goal)
        newE = next_ecut(cfg, Ecut)
        bres = newE > Ecut + 1e-10 ? basis_estimate!(refcache, used_ref, spent, cfg, current, newE, goal) : (; estimate=0.0, wall=0.0)
        sest = scf_anchor * (tol / max(cfg.start_tol, eps()))
        total_est = bres.estimate + k_est + sest
        basis_cost = newE > Ecut + 1e-10 ? bres.wall + current.wall * max((newE / Ecut)^1.5 - 1, 0.05) : Inf
        nextnk = min(max_nk, nk + 1)
        k_cost = nk < max_nk ? (k_candidate !== nothing && k_candidate_n == nextnk ? k_candidate.wall : current.wall * max((nextnk / nk)^3, 1.15)) : Inf
        scf_cost = tol > 1.01e-9 ? current.wall : Inf
        stop_threshold = SAFETY_FACTOR * target
        action = total_est <= stop_threshold ? :stop : choose_action(bres.estimate, k_est, sest, basis_cost, k_cost, scf_cost, Ecut, newE, nk, max_nk, tol)
        verified_stop = false
        if action == :kgrid && total_est <= NEAR_THRESHOLD_FACTOR * stop_threshold && nk < max_nk && k_source != "direct_probe" && k_source != "near_threshold_verification"
            nextnk = nk + 1
            probe = get_run!(run_cache, used, spent, cfg, Ecut, nextnk, tol, goal)
            verified = NEAR_THRESHOLD_VERIFY_INFLATION * property_difference(current, probe, cfg, goal)
            k_est = verified
            k_source = "near_threshold_verification"
            k_candidate = probe
            k_candidate_n = nextnk
            total_est = bres.estimate + k_est + sest
            if total_est <= stop_threshold
                action = :stop
                verified_stop = true
            else
                action = choose_action(bres.estimate, k_est, sest, basis_cost, probe.wall, scf_cost, Ecut, newE, nk, max_nk, tol)
            end
        end
        if action == :stop && nk <= cfg.ref_nk && k_source != "direct_probe" && k_source != "stop_verification" && k_source != "near_threshold_verification"
            nextnk = nk + 1
            probe = get_run!(run_cache, used, spent, cfg, Ecut, nextnk, tol, goal)
            verified = property_difference(current, probe, cfg, goal)
            k_est = verified
            k_source = "stop_verification"
            k_candidate = probe
            k_candidate_n = nextnk
            total_est = bres.estimate + k_est + sest
            verified_stop = true
            if total_est > stop_threshold
                action = choose_action(bres.estimate, k_est, sest, basis_cost, probe.wall, scf_cost, Ecut, newE, nk, max_nk, tol)
            end
        end
        push!(trace, (cfg.name, string(goal), target, step, Ecut, nk, tol, bres.estimate, k_est, sest, total_est, actual, string(action), spent[], k_source, "anchored_probe", verified_stop))
        action == :stop && break
        if action == :Ecut
            Ecut = newE
            current = get_run!(run_cache, used, spent, cfg, Ecut, nk, tol, goal)
            k_est = isfinite(k_est) ? 1.10 * k_est : k_est
            k_source = "transferred_prediction"
            k_candidate = nothing
        elseif action == :kgrid
            old_nk = nk
            nextnk = nk + 1
            if k_candidate === nothing || k_candidate_n != nextnk || abs(k_candidate.Ecut - Ecut) > 1e-10 || abs(k_candidate.tol - tol) > 1e-15
                k_candidate = get_run!(run_cache, used, spent, cfg, Ecut, nextnk, tol, goal)
                k_candidate_n = nextnk
            end
            observed = property_difference(current, k_candidate, cfg, goal)
            current = k_candidate
            nk = nextnk
            k_est = predict_next_k_error(observed, old_nk, nk)
            k_source = nk < cfg.ref_nk ? "decay_prediction" : (nk == cfg.ref_nk ? "boundary_prediction" : "boundary_tail_prediction")
            k_candidate = nothing
        elseif action == :scf
            tol2 = next_tol(tol)
            tight = get_run!(run_cache, used, spent, cfg, Ecut, nk, tol2, goal)
            observed = property_difference(current, tight, cfg, goal)
            current = tight
            tol = tol2
            scf_anchor = observed
            k_candidate = nothing
        else
            break
        end
    end
    actual = property_difference(current, truth, cfg, goal)
    success = actual <= target
    outcome = (cfg.name, string(goal), target, success, actual, spent[], current.Ecut, current.nk, current.tol, length(trace), stringent_wall)
    trace, outcome
end

function uniform_goal(cfg, target, goal, truth, run_cache)
    used = Set{Any}()
    spent = Ref(0.0)
    Ecut = cfg.start_Ecut
    nk = cfg.start_nk
    tol = cfg.start_tol
    current = get_run!(run_cache, used, spent, cfg, Ecut, nk, tol, goal)
    rows = Any[]
    for step in 1:MAX_CONTROL_STEPS
        newE = next_ecut(cfg, Ecut)
        newnk = min(cfg.ref_nk, nk + 1)
        newtol = next_tol(tol)
        if newE == Ecut && newnk == nk && newtol == tol
            break
        end
        nxt = get_run!(run_cache, used, spent, cfg, newE, newnk, newtol, goal)
        successive = property_difference(current, nxt, cfg, goal)
        actual = property_difference(nxt, truth, cfg, goal)
        push!(rows, (cfg.name, string(goal), target, step, newE, newnk, newtol, successive, actual, spent[]))
        current = nxt
        Ecut, nk, tol = newE, newnk, newtol
        successive <= SAFETY_FACTOR * target && break
    end
    actual = property_difference(current, truth, cfg, goal)
    outcome = (cfg.name, string(goal), target, actual <= target, actual, spent[], current.Ecut, current.nk, current.tol, length(rows))
    rows, outcome
end

function warmup()
    cfg = (name="warmup", geom=diamond(:Si, 5.431; displacement=[0.005,0.0,0.0]), start_Ecut=10.0, ref_Ecut=13.0, truth_Ecut=15.0, start_nk=1, ref_nk=1, truth_nk=1, start_tol=1e-3)
    try
        r = run_scf(cfg, 10.0, 1, 1e-2, :force)
        b = PlaneWaveBasis(r.scfres.basis.model; Ecut=13.0, kgrid=[1,1,1], architecture=DFTK.CPU())
        ref = refine_scfres(r.scfres, b; tol=1e-3, callback=identity)
        refine_energies(ref)
        refine_forces(ref)
    catch
    end
    GC.gc()
end

function main()
    println("Precompiling DFTK, force and refinement kernels; warmup is excluded from reported timings.")
    warmup()
    write_log_header()
    tmain = time_ns()
    append_log("event", "start", nowstamp())
    energy_trace = Any[]
    energy_outcomes = Any[]
    uniform_rows = Any[]
    uniform_outcomes = Any[]
    system_rows = Any[]
    reference_rows = Any[]
    for cfg in energy_systems()
        append_log("reference_start", cfg.name, "energy", nowstamp())
        ref = reference_run(cfg, :energy)
        truth = truth_run(cfg, :energy)
        referr = property_difference(ref, truth, cfg, :energy)
        push!(system_rows, (cfg.name, "energy", length(cfg.geom[2]), cfg.start_Ecut, cfg.ref_Ecut, cfg.truth_Ecut, cfg.start_nk, cfg.ref_nk, cfg.truth_nk, cfg.start_tol, ref.wall, truth.wall, ref.energy, truth.energy))
        push!(reference_rows, (cfg.name, "energy", referr, ref.wall, truth.wall, cfg.ref_Ecut, cfg.truth_Ecut, cfg.ref_nk, cfg.truth_nk))
        append_log("reference_validation", cfg.name, "energy", "difference", referr)
        run_cache = Dict{Any,Any}()
        refcache = Dict{Any,Any}()
        for target in ENERGY_TARGETS_MEV
            tr, out = adaptive_goal(cfg, target, :energy, truth, ref.wall, run_cache, refcache)
            append!(energy_trace, tr)
            push!(energy_outcomes, out)
            ur, uo = uniform_goal(cfg, target, :energy, truth, run_cache)
            append!(uniform_rows, ur)
            push!(uniform_outcomes, uo)
            append_log("target_done", cfg.name, "energy", target, "success", out[4], "actual_error", out[5], "adaptive_cost_s", out[6])
        end
    end
    force_trace = Any[]
    force_outcomes = Any[]
    force_uniform_rows = Any[]
    force_uniform_outcomes = Any[]
    for cfg in force_systems()
        append_log("reference_start", cfg.name, "force", nowstamp())
        ref = reference_run(cfg, :force)
        truth = truth_run(cfg, :force)
        referr = property_difference(ref, truth, cfg, :force)
        push!(system_rows, (cfg.name, "force", length(cfg.geom[2]), cfg.start_Ecut, cfg.ref_Ecut, cfg.truth_Ecut, cfg.start_nk, cfg.ref_nk, cfg.truth_nk, cfg.start_tol, ref.wall, truth.wall, ref.energy, truth.energy))
        push!(reference_rows, (cfg.name, "force", referr, ref.wall, truth.wall, cfg.ref_Ecut, cfg.truth_Ecut, cfg.ref_nk, cfg.truth_nk))
        append_log("reference_validation", cfg.name, "force", "difference", referr)
        run_cache = Dict{Any,Any}()
        refcache = Dict{Any,Any}()
        for target in FORCE_TARGETS_EVANG
            tr, out = adaptive_goal(cfg, target, :force, truth, ref.wall, run_cache, refcache)
            append!(force_trace, tr)
            push!(force_outcomes, out)
            ur, uo = uniform_goal(cfg, target, :force, truth, run_cache)
            append!(force_uniform_rows, ur)
            push!(force_uniform_outcomes, uo)
            append_log("target_done", cfg.name, "force", target, "success", out[4], "actual_error", out[5], "adaptive_cost_s", out[6])
        end
    end
    force_compare_rows = Any[]
    for out in force_outcomes
        system, goal, target = out[1], out[2], out[3]
        idx = findfirst(u -> u[1] == system && u[2] == goal && abs(u[3] - target) < 1e-12, force_uniform_outcomes)
        idx === nothing && error("Missing force uniform comparison for $(system), target=$(target)")
        uo = force_uniform_outcomes[idx]
        speedup = out[6] > 0 ? uo[6] / out[6] : NaN
        push!(force_compare_rows, (out..., uo[4], uo[5], uo[6], speedup))
    end
    write_tsv(joinpath(MAIN_DIR, "fig1_estimator_calibration.dat"),
        ["system","goal","target","step","Ecut_Ha","kgrid_n","scf_tol","basis_estimate","kpoint_estimate","scf_estimate","estimated_controllable_error","actual_error","selected_action","cumulative_cost_s","kpoint_estimate_source","scf_estimate_source","stop_verified"], energy_trace)
    energy_compare_rows = Any[]
    for out in energy_outcomes
        system, goal, target = out[1], out[2], out[3]
        idx = findfirst(u -> u[1] == system && u[2] == goal && abs(u[3] - target) < 1e-12, uniform_outcomes)
        idx === nothing && error("Missing uniform comparison for $(system), target=$(target)")
        uo = uniform_outcomes[idx]
        speedup = out[6] > 0 ? uo[6] / out[6] : NaN
        push!(energy_compare_rows, (out..., uo[4], uo[5], uo[6], speedup))
    end
    write_tsv(joinpath(MAIN_DIR, "fig2_cost_to_target.dat"),
        ["system","goal","target_meV_per_atom","success","actual_error_meV_per_atom","adaptive_cost_s","final_Ecut_Ha","final_kgrid_n","final_scf_tol","control_steps","stringent_reference_time_s","uniform_success","uniform_actual_error_meV_per_atom","uniform_cost_s","uniform_over_adaptive_speedup"], energy_compare_rows)
    write_tsv(joinpath(MAIN_DIR, "fig3_error_budget.dat"),
        ["system","goal","target","step","Ecut_Ha","kgrid_n","scf_tol","basis_estimate","kpoint_estimate","scf_estimate","estimated_controllable_error","actual_error","selected_action","cumulative_cost_s","kpoint_estimate_source","scf_estimate_source","stop_verified"], energy_trace)
    write_tsv(joinpath(MAIN_DIR, "fig4_parameter_path.dat"),
        ["system","goal","target","step","Ecut_Ha","kgrid_n","scf_tol","basis_estimate","kpoint_estimate","scf_estimate","estimated_controllable_error","actual_error","selected_action","cumulative_cost_s","kpoint_estimate_source","scf_estimate_source","stop_verified"], energy_trace)
    write_tsv(joinpath(MAIN_DIR, "fig5_target_reliability.dat"),
        ["system","goal","target_meV_per_atom","success","actual_error_meV_per_atom","adaptive_cost_s","final_Ecut_Ha","final_kgrid_n","final_scf_tol","control_steps","stringent_reference_time_s"], energy_outcomes)
    write_tsv(joinpath(MAIN_DIR, "table1_systems.dat"),
        ["system","goal","atoms","start_Ecut_Ha","reference_Ecut_Ha","truth_Ecut_Ha","start_kgrid_n","reference_kgrid_n","truth_kgrid_n","start_scf_tol","reference_wall_time_s","truth_wall_time_s","reference_energy_Ha","truth_energy_Ha"], system_rows)
    write_tsv(joinpath(MAIN_DIR, "table2_energy_targets.dat"),
        ["system","goal","target_meV_per_atom","success","actual_error_meV_per_atom","adaptive_cost_s","final_Ecut_Ha","final_kgrid_n","final_scf_tol","control_steps","stringent_reference_time_s","uniform_success","uniform_actual_error_meV_per_atom","uniform_cost_s","uniform_over_adaptive_speedup"], energy_compare_rows)
    write_tsv(joinpath(SI_DIR, "figS1_uniform_protocol.dat"),
        ["system","goal","target","step","Ecut_Ha","kgrid_n","scf_tol","successive_change","actual_error","cumulative_cost_s"], uniform_rows)
    write_tsv(joinpath(SI_DIR, "tableS1_uniform_protocol.dat"),
        ["system","goal","target","success","actual_error","uniform_cost_s","final_Ecut_Ha","final_kgrid_n","final_scf_tol","steps"], uniform_outcomes)
    write_tsv(joinpath(SI_DIR, "figS2_force_estimator_calibration.dat"),
        ["system","goal","target","step","Ecut_Ha","kgrid_n","scf_tol","basis_estimate","kpoint_estimate","scf_estimate","estimated_controllable_error","actual_error","selected_action","cumulative_cost_s","kpoint_estimate_source","scf_estimate_source","stop_verified"], force_trace)
    write_tsv(joinpath(SI_DIR, "figS3_force_cost.dat"),
        ["system","goal","target_eV_per_A","success","actual_error_eV_per_A","adaptive_cost_s","final_Ecut_Ha","final_kgrid_n","final_scf_tol","control_steps","stringent_reference_time_s","uniform_success","uniform_actual_error_eV_per_A","uniform_cost_s","uniform_over_adaptive_speedup"], force_compare_rows)
    write_tsv(joinpath(SI_DIR, "figS4_force_uniform.dat"),
        ["system","goal","target","step","Ecut_Ha","kgrid_n","scf_tol","successive_change","actual_error","cumulative_cost_s"], force_uniform_rows)
    write_tsv(joinpath(SI_DIR, "figS5_energy_trajectories.dat"),
        ["system","goal","target","step","Ecut_Ha","kgrid_n","scf_tol","basis_estimate","kpoint_estimate","scf_estimate","estimated_controllable_error","actual_error","selected_action","cumulative_cost_s","kpoint_estimate_source","scf_estimate_source","stop_verified"], energy_trace)
    write_tsv(joinpath(SI_DIR, "figS6_terminal_budget.dat"),
        ["system","goal","target","step","Ecut_Ha","kgrid_n","scf_tol","basis_estimate","kpoint_estimate","scf_estimate","estimated_controllable_error","actual_error","selected_action","cumulative_cost_s","kpoint_estimate_source","scf_estimate_source","stop_verified"], energy_trace)
    write_tsv(joinpath(SI_DIR, "figS7_final_parameters.dat"),
        ["system","goal","target_meV_per_atom","success","actual_error_meV_per_atom","adaptive_cost_s","final_Ecut_Ha","final_kgrid_n","final_scf_tol","control_steps","stringent_reference_time_s"], energy_outcomes)
    write_tsv(joinpath(SI_DIR, "figS8_reference_validation.dat"),
        ["system","goal","reference_vs_truth_error","reference_wall_time_s","truth_wall_time_s","reference_Ecut_Ha","truth_Ecut_Ha","reference_kgrid_n","truth_kgrid_n"], reference_rows)
    write_tsv(joinpath(SI_DIR, "tableS2_force_targets.dat"),
        ["system","goal","target_eV_per_A","success","actual_error_eV_per_A","adaptive_cost_s","final_Ecut_Ha","final_kgrid_n","final_scf_tol","control_steps","stringent_reference_time_s","uniform_success","uniform_actual_error_eV_per_A","uniform_cost_s","uniform_over_adaptive_speedup"], force_compare_rows)
    write_tsv(joinpath(SI_DIR, "tableS3_force_uniform.dat"),
        ["system","goal","target","success","actual_error","uniform_cost_s","final_Ecut_Ha","final_kgrid_n","final_scf_tol","steps"], force_uniform_outcomes)
    write_tsv(joinpath(SI_DIR, "tableS4_reference_validation.dat"),
        ["system","goal","reference_vs_truth_error","reference_wall_time_s","truth_wall_time_s","reference_Ecut_Ha","truth_Ecut_Ha","reference_kgrid_n","truth_kgrid_n"], reference_rows)
    append_log("event", "finish", nowstamp(), "production_wall_s", (time_ns()-tmain)/1e9)
    println("GoalDFT.jl completed.")
    println("Run: python redraw_paper.py")
end

main()

