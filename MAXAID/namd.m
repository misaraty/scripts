clc;
clear;
tic;
[filepath] = fileparts([mfilename('fullpath'), '.m']);
addpath(genpath(filepath));
cd(filepath);
format short e;
global active_space hbar kb Ry_to_eV Ha_to_eV rt nucl_dt elec_dt integrator sh_algo num_sh_traj decoherence_method boltz_flag Temp active_map outdir;
rt = [fullfile(filepath, '..'), filesep];
outdir = fullfile(filepath, 'out');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
namdtime = 10000;
num_sh_traj = 1000;
iconds_i = 1000;
Ham_size = 4001;
active_space = [1, 2];
states = {{'GS', [1, -1], 0.00}, {'S1', [1, -2], 0.00}};
decoherence_method = 2;
nucl_dt = 1;
elec_dt = 1;
integrator = 0;
sh_algo = 0;
boltz_flag = 1;
Temp = 300;
alp_bet = 0;
hbar = 0.658218;
kb = 8.617e-5;
Ry_to_eV = 13.60569253;
Ha_to_eV = 27.211396;
me_numstates = numel(states);
numstates = numel(active_space);
me_states = cell(1, me_numstates);
for me_state_i = 1:me_numstates
    me_states{me_state_i} = me_state(states{me_state_i});
end
num_elec = numel(me_states{1}.actual_state);
active_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
for k = 1:numstates
    active_map(active_space(k)) = k - 1;
end
iconds = zeros(iconds_i, 2);
for ii = 1:iconds_i
    iconds(ii, :) = [ii - 1, me_numstates - 1];
end
H_batch = cell(Ham_size, 1);
for j = 0:(Ham_size - 1)
    Ham_re = importdata(fullfile(rt, 'res', ['0_Ham_', num2str(j), '_re']));
    Ham_im = importdata(fullfile(rt, 'res', ['0_Ham_', num2str(j), '_im']));
    Ham = extract_2D(Ham_re, active_space, -1) + 1i * extract_2D(Ham_im, active_space, -1);
    H_batch{j + 1} = Ham * Ry_to_eV;
end
for icond = 0:(size(iconds, 1) - 1)
    oe_es = cell(namdtime, 1);
    me_es = cell(namdtime, 1);
    for celln = 1:namdtime
        oe_es{celln} = ElectronicStructure(2*numstates);
        me_es{celln} = ElectronicStructure(me_numstates);
    end
    istart = iconds(icond+1, 1);
    for j = istart:(istart + namdtime - 1)
        t = j - istart + 1;
        j1 = mod(j, Ham_size);
        Hij = H_batch{j1+1};
        Htmp = zeros(2*numstates, 2*numstates);
        for k1 = 0:(numstates - 1)
            i1 = 2 * k1 + 1;
            i2 = 2 * k1 + 2;
            for k2 = 0:(numstates - 1)
                jx1 = 2 * k2 + 1;
                jx2 = 2 * k2 + 2;
                Htmp(i1, jx1) = Hij(k1+1, k2+1);
                Htmp(i2, jx2) = Hij(k1+1, k2+1);
            end
        end
        if alp_bet ~= 0
            for k1 = 0:(numstates - 1)
                i1 = 2 * k1 + 1;
                i2 = 2 * k1 + 2;
                for k2 = 0:(numstates - 1)
                    jx1 = 2 * k2 + 1;
                    jx2 = 2 * k2 + 2;
                    Htmp(i2, jx1) = Hij(k1+1, k2+1);
                    Htmp(i1, jx2) = Hij(k1+1, k2+1);
                end
            end
        end
        oe_es{t}.Hcurr = Htmp;
    end
    me_es{1} = ElectronicStructure_set_state(me_es{1}, iconds(icond+1, 2));
    for j = istart:(istart + namdtime - 1)
        t = j - istart + 1;
        Hme = zeros(me_numstates, me_numstates);
        for I = 0:(me_numstates - 1)
            for J = 0:(me_numstates - 1)
                orb_i = 0;
                orb_j = 0;
                [delt, ~, ~, orb_i, orb_j] = delta(me_states{I+1}.actual_state, me_states{J+1}.actual_state, orb_i, orb_j);
                if delt
                    orb_i = ext2int(orb_i, active_space);
                    orb_j = ext2int(orb_j, active_space);
                    Hme(I+1, J+1) = Hme(I+1, J+1) + oe_es{t}.Hcurr(orb_i+1, orb_j+1);
                end
            end
            diagE = 0.0;
            for el = 0:(num_elec - 1)
                orb_i = ext2int(me_states{I+1}.actual_state(el+1), active_space);
                diagE = diagE + oe_es{t}.Hcurr(orb_i+1, orb_i+1);
            end
            shift = me_states{I+1}.shift;
            Hme(I+1, I+1) = diagE + shift;
        end
        me_es{t}.Hcurr = Hme;
    end
    outfile = fullfile(outdir, ['me_energies', num2str(icond)]);
    out = fopen(outfile, 'w');
    for j = istart:(istart + namdtime - 1)
        t = j - istart + 1;
        fprintf(out, 't= %d  E[0]= %.12g  ', j, real(me_es{t}.Hcurr(1, 1)));
        e0 = me_es{t}.Hcurr(1, 1);
        for I = 0:(me_es{t}.num_states - 1)
            fprintf(out, 'E[%d]-E[0]= %.12g  ', I, real(me_es{t}.Hcurr(I+1, I+1)-e0));
        end
        fprintf(out, '\n');
    end
    fclose(out);
    if decoherence_method > 0
        me_es = run_decoherence_rates(me_es, icond);
    end
    [me_es, me_states] = run_namd(me_es, me_states, icond);
    oe_es = [];
    me_es = [];
end
elapsedTime = toc;
logFile = fopen(fullfile(filepath, 'log.txt'), 'w');
fprintf('Running time: %.1f seconds\n', elapsedTime);
fprintf(logFile, 'Running time: %.1f seconds\n', elapsedTime);
fclose(logFile);
fileID = fopen(fullfile(outdir, 'me_pop0'), 'r');
last_line = '';
while ~feof(fileID)
    current_line = fgetl(fileID);
    if ischar(current_line)
        last_line = current_line;
    end
end
fclose(fileID);
disp(last_line);
function this = ElectronicStructure(n)
this.curr_state = 0;
this.num_states = n;
this.Ccurr = zeros(n, 1);
this.Cprev = zeros(n, 1);
this.Cnext = zeros(n, 1);
this.Hcurr = zeros(n, n);
this.Hprev = zeros(n, n);
this.Hnext = zeros(n, n);
this.dHdt = zeros(n, n);
this.g = eye(n);
this.tau_m = zeros(n, 1);
this.t_m = zeros(n, 1);
this.A = zeros(n, n);
end

function out = extract_2D(in, templ, shift)
idx = templ + shift + 1;
out = in(idx, idx);
end

function this = me_state(cs_)
this.name = cs_{1};
this.actual_state = cs_{2};
if length(cs_) >= 3
    this.shift = cs_{3};
else
    this.shift = 0.0;
end
end

function [res, A, B, a, b] = delta(A, B, a, b)
C_ = unique([A, B]);
nexc = 0;
for i = 1:numel(C_)
    n_in_a = sum(A == C_(i));
    n_in_b = sum(B == C_(i));
    d = n_in_a - n_in_b;
    if d > 0
        nexc = nexc + d;
    end
    if d == 1
        a = C_(i);
    end
    if d == -1
        b = C_(i);
    end
end
res = (nexc == 1);
end

function internal = ext2int(external, active_space)
global active_map;
idx = active_map(abs(external));
f = double(external < 0);
internal = 2 * abs(idx) + f;
end

function [X, Y, a, b] = regression(X, Y, a, b)
sx = sum(X);
if abs(sx) < eps
    b = 0.0;
else
    b = sum(Y) / sx;
end
end

function [result, x] = decoherence_rates(x, dt, icond, i1, j1)
global hbar outdir;
len = numel(x);
sz = floor(len/2);
C = zeros(sz, 1);
IC = zeros(sz, 1);
IIC = zeros(sz, 1);
D = zeros(sz, 1);
T = zeros(sz, 1);
selIIC = zeros(sz, 1);
for t = 0:(sz - 1)
    C(t+1) = sum(x(1:sz).*x((1 + t):(sz + t))) / sz;
end
sum0 = 0.0;
for t = 0:(sz - 1)
    IC(t+1) = sum0;
    sum0 = sum0 + C(t+1) * (dt / hbar);
end
sum0 = 0.0;
for t = 0:(sz - 1)
    IIC(t+1) = sum0;
    sum0 = sum0 + IC(t+1) * (dt / hbar);
end
D = exp(-IIC);
nrm = C(1);
if abs(nrm) > eps
    C = C / nrm;
end
dE = 0.0025;
Npoints = 2000;
J = zeros(Npoints, 1);
tt = (1:(sz - 1))' * dt;
Ct = C(2:end);
out1 = fopen(fullfile(outdir, ['icond', num2str(icond), 'pair', num2str(i1), '_', num2str(j1), 'Spectral_density.txt']), 'w');
for w = 0:(Npoints - 1)
    ww = w * dE;
    J(w+1) = (1.0 + 2.0 * sum(cos(ww*tt).*Ct)) * dt;
    J(w+1) = J(w+1) * J(w+1) / (2.0 * pi);
    fprintf(out1, 'w(eV)= %.12g w(cm^-1)= %.12g J= %.12g sqrt(J)= %.12g\n', w*dE, w*dE*8065.54468111324, J(w+1), sqrt(J(w+1)));
end
fclose(out1);
first = true;
cnt = 0;
for t = 0:(sz - 1)
    if first
        if IIC(t+1) < 2.3
            cnt = cnt + 1;
            T(cnt) = t * t * dt * dt;
            selIIC(cnt) = IIC(t+1);
        else
            first = false;
        end
    end
end
T = T(1:cnt);
selIIC = selIIC(1:cnt);
a = 0;
b = 0;
if ~isempty(T)
    [~, ~, a, b] = regression(T, selIIC, a, b);
end
if b < 0.0
    b = 0.0;
end
out = fopen(fullfile(outdir, ['icond', num2str(icond), 'pair', num2str(i1), '_', num2str(j1), 'Dephasing_function.txt']), 'w');
fprintf(out, 'Time    D(t)       fitted D(t)     Normalized_autocorrelation_function  Unnormalized_autocorrelation_function   Second cumulant\n');
for t = 0:(sz - 1)
    fprintf(out, '%.12g  %.12g  %.12g  %.12g %.12g  %.12g\n', t*dt, D(t+1), exp(-a)*exp(-b*t*t*dt*dt), C(t+1), nrm*C(t+1), IIC(t+1));
end
fclose(out);
result = sqrt(b);
end

function me_es = run_decoherence_rates(me_es, icond)
global nucl_dt outdir;
sz = numel(me_es);
N = me_es{1}.num_states;
rij = zeros(N, N);
out = fopen(fullfile(outdir, ['decoherence_rates_icond', num2str(icond), '.txt']), 'w');
for i1 = 0:(N - 1)
    for j1 = 0:(N - 1)
        if i1 == j1
            rij(i1+1, j1+1) = 0.0;
        else
            Eij = zeros(sz, 1);
            ave_dEij = 0.0;
            for t = 0:(sz - 1)
                dEij = real(me_es{t+1}.Hcurr(i1+1, i1+1)-me_es{t+1}.Hcurr(j1+1, j1+1));
                Eij(t+1) = dEij;
                ave_dEij = ave_dEij + dEij;
            end
            Eij = Eij - ave_dEij / sz;
            [rij(i1+1, j1+1), ~] = decoherence_rates(Eij, nucl_dt, icond, i1, j1);
        end
        fprintf(out, '%.12g ', real(rij(i1+1, j1+1)));
    end
    fprintf(out, '\n');
end
fclose(out);
end

function this = ElectronicStructure_set_state(this, indx)
this.Ccurr(:) = 0;
this.Ccurr(indx+1) = 1 + 0i;
this.curr_state = indx;
end

function this = ElectronicStructure_insertion(this, es)
this.num_states = es.num_states;
this.curr_state = es.curr_state;
this.Cprev = es.Cprev;
this.Ccurr = es.Ccurr;
this.Cnext = es.Cnext;
this.tau_m = es.tau_m;
this.t_m = es.t_m;
this.A = es.A;
end

function this = ElectronicStructure_init_hop_prob1(this)
this.g = eye(this.num_states);
end

function [E, Eex] = Efield(t, E, Eex)
Eex = 0.0;
E = 0.0;
end

function this = ElectronicStructure_rot1(this, phi, i, j)
c = cos(phi);
s = sin(phi);
ci = this.Ccurr(i + 1);
cj = this.Ccurr(j + 1);
this.Ccurr(i + 1) = c * ci + s * cj;
this.Ccurr(j + 1) = -s * ci + c * cj;
end

function this = ElectronicStructure_rot2(this, phi, i, j)
cs = cos(phi);
isi = 1i * sin(phi);
ci = this.Ccurr(i + 1);
cj = this.Ccurr(j + 1);
this.Ccurr(i + 1) = cs * ci + isi * cj;
this.Ccurr(j + 1) = isi * ci + cs * cj;
end

function this = ElectronicStructure_rot(this, Hij, dt, i, j)
global hbar;
phi1 = 0.5 * dt * imag(Hij) / hbar;
phi2 = -dt * real(Hij) / hbar;
this = ElectronicStructure_rot1(this, phi1, i, j);
this = ElectronicStructure_rot2(this, phi2, i, j);
this = ElectronicStructure_rot1(this, phi1, i, j);
end

function this = ElectronicStructure_phase(this, Hii, dt, i)
global hbar;
phi = -dt * real(Hii) / hbar;
this.Ccurr(i + 1) = exp(1i*phi) * this.Ccurr(i + 1);
end

function this = ElectronicStructure_propagate_coefficients(this, dt, Ef)
for i = 0:(this.num_states - 1)
    for j = (i + 1):(this.num_states - 1)
        this = ElectronicStructure_rot(this, this.Hcurr(i + 1, j + 1), 0.5*dt, i, j);
    end
end
for i = 0:(this.num_states - 1)
    this = ElectronicStructure_phase(this, this.Hcurr(i + 1, i + 1), dt, i);
end
for i = (this.num_states - 1):-1:0
    for j = (this.num_states - 1):-1:(i + 1)
        this = ElectronicStructure_rot(this, this.Hcurr(i + 1, j + 1), 0.5*dt, i, j);
    end
end
end

function this = ElectronicStructure_update_populations(this)
this.A = conj(this.Ccurr) * (this.Ccurr.');
end

function this = ElectronicStructure_update_hop_prob_fssh(this, dt, boltz_flag, Temp, Ef, Eex, rates)
global hbar kb;
this = ElectronicStructure_update_populations(this);
Heff = this.Hcurr;
for i = 0:(this.num_states - 1)
    a_ii = real(this.A(i + 1, i + 1));
    if a_ii < 1e-12
        a_ii = 1e-12;
    end
    g_row_sum = 0.0;
    for j = 0:(this.num_states - 1)
        if j ~= i
            g_ij = (2.0 * dt / (a_ii * hbar)) * imag(this.A(i + 1, j + 1)*Heff(i + 1, j + 1));
            if g_ij < 0.0
                g_ij = 0.0;
            end
            E_i = real(Heff(i + 1, i + 1));
            E_j = real(Heff(j + 1, j + 1));
            dE = E_j - E_i;
            bf = 1.0;
            if dE > Eex
                bf = exp(-((dE - Eex) / (kb * Temp)));
            end
            g_ij = g_ij * bf;
            g_row_sum = g_row_sum + g_ij;
            this.g(i + 1, j + 1) = g_ij;
        end
    end
    this.g(i + 1, i + 1) = 1.0 - g_row_sum;
end
end

function [es, rates] = propagate_electronic(es, i, rates)
global integrator sh_algo nucl_dt elec_dt boltz_flag Temp;
nel = nucl_dt / elec_dt;
tim = 0;
Eex = 0;
Ef = zeros(3, 1);
if integrator == 0
    for j = 0:(nel - 1)
        tim = i * nucl_dt + j * elec_dt;
        [Ef, Eex] = Efield(tim, Ef, Eex);
        es{i + 1} = ElectronicStructure_propagate_coefficients(es{i + 1}, elec_dt, Ef);
        if sh_algo == 0
            es{i + 1} = ElectronicStructure_update_hop_prob_fssh(es{i + 1}, elec_dt, boltz_flag, Temp, Ef, Eex, rates);
        end
    end
end
end

function [sh_prob, hopstate] = hop(sh_prob, hopstate, numstates)
in = hopstate;
ksi = rand();
probs = sh_prob(in+1, :);
nrm = sum(probs);
if nrm <= 0
    return;
end
probs = probs / nrm;
cprobs = cumsum(probs);
hstate = find((0 < cprobs) & (ksi <= cprobs), 1, 'first') - 1;
if isempty(hstate)
    hstate = in;
end
hopstate = hstate;
end

function this = ElectronicStructure_update_decoherence_times(this, rates)
this = ElectronicStructure_update_populations(this);
pops = real(diag(this.A));
this.tau_m = real(rates) * pops;
end

function this = ElectronicStructure_decohere(this, i)
this.Ccurr(:) = 0;
this.Ccurr(i + 1) = 1.0;
this.curr_state = i;
end

function this = ElectronicStructure_project_out(this, i)
this = ElectronicStructure_update_populations(this);
this.Ccurr(i + 1) = 0.0;
nrm2 = sum(real(diag(this.A))) - real(this.A(i + 1, i + 1));
if nrm2 <= 0
    this.Ccurr(:) = 0.0;
    this.Ccurr(this.curr_state+1) = 1.0;
else
    this.Ccurr = this.Ccurr / sqrt(nrm2);
end
this = ElectronicStructure_update_populations(this);
end

function this = ElectronicStructure_dish_decoherence(this, dt, boltz_flag, Temp, rates)
global kb;
this = ElectronicStructure_update_decoherence_times(this, rates);
for i = 0:(this.num_states - 1)
    if this.tau_m(i + 1) <= 0
        continue;
    end
    rnd_i = 1.0 / this.tau_m(i + 1);
    if this.t_m(i + 1) >= rnd_i
        zeta = rand();
        P = real(this.A(i + 1, i + 1));
        dE = real(this.Hcurr(i + 1, i + 1)-this.Hcurr(this.curr_state+1, this.curr_state+1));
        if dE > 0
            P = P * exp(-(dE / (kb * Temp)));
        end
        if zeta < P
            this = ElectronicStructure_decohere(this, i);
            break;
        else
            this = ElectronicStructure_project_out(this, i);
            this.t_m(i + 1) = 0;
            this.tau_m(i + 1) = 0.0;
        end
    end
end
this.t_m = this.t_m + dt;
end

function this = ElectronicStructure_sdm_decoherence(this, dt, act_st, rates, tol)
if nargin < 5
    tol = 0.0;
end
if act_st < 0 || act_st >= this.num_states
    error('Error in ElectronicStructure_sdm_decoherence: active state index out of range');
end
this = ElectronicStructure_update_populations(this);
p_aa_old = real(this.A(act_st+1, act_st+1));
if p_aa_old > 1.0 + tol
    sclf = 1.0 / sqrt(p_aa_old);
    this.Ccurr = this.Ccurr * sclf;
    this = ElectronicStructure_update_populations(this);
    p_aa_old = real(this.A(act_st+1, act_st+1));
end
if p_aa_old <= 0.0
    return;
end
inact_st_pop = 0.0;
for i = 0:(this.num_states - 1)
    if i ~= act_st
        itau = real(rates(i + 1, act_st + 1));
        sclf = exp(-dt*itau);
        this.Ccurr(i + 1) = this.Ccurr(i + 1) * sclf;
        inact_st_pop = inact_st_pop + abs(this.Ccurr(i + 1))^2;
    end
end
if inact_st_pop > 1.0
    error('Error in ElectronicStructure_sdm_decoherence: inactive-state population > 1.0');
end
p_aa_new = 1.0 - inact_st_pop;
if p_aa_new < 0.0
    error('Error in ElectronicStructure_sdm_decoherence: new active-state population is negative');
end
sclf = sqrt(p_aa_new/p_aa_old);
this.Ccurr(act_st+1) = this.Ccurr(act_st+1) * sclf;
this = ElectronicStructure_update_populations(this);
new_norm = real(this.Ccurr'*this.Ccurr);
if abs(new_norm-1.0) > 0.1
    error('Error in ElectronicStructure_sdm_decoherence: norm deviates too much from 1');
end
end

function this = ElectronicStructure_gfsh_decohere(this, dt, Temp)
global kb;
this = ElectronicStructure_update_populations(this);
curr_state = this.curr_state;
pop_diag = real(diag(this.A));
total_pop = sum(pop_diag);
if total_pop <= 0
    return;
end
normalized_pop = pop_diag / total_pop;
num_decohere = sum(normalized_pop > rand(size(normalized_pop)));
if num_decohere <= 0
    this.t_m = this.t_m + dt;
    return;
end
decohere_states = randperm(this.num_states, min(num_decohere, this.num_states)) - 1;
for idx = 1:numel(decohere_states)
    i = decohere_states(idx);
    zeta = rand();
    P = real(this.A(i + 1, i + 1)) / total_pop;
    dE = real(this.Hcurr(i + 1, i + 1)-this.Hcurr(curr_state+1, curr_state+1));
    if dE > 0
        P = P * exp(-(dE / (kb * Temp)));
    end
    if zeta < P
        this = ElectronicStructure_decohere(this, i);
        curr_state = i;
    else
        this = ElectronicStructure_project_out(this, i);
    end
    this.t_m(i + 1) = 0;
    this.tau_m(i + 1) = 0.0;
end
this.t_m = this.t_m + dt;
end

function [me_es, me_states] = run_namd(me_es, me_states, icond)
global nucl_dt elec_dt num_sh_traj decoherence_method outdir Temp;
sz = numel(me_es);
nst = me_es{1}.num_states;
init_state = me_es{1}.curr_state;
sh_pops = zeros(sz, nst);
se_pops = zeros(sz, nst);
rates = zeros(nst, nst);
if decoherence_method > 0
    filename = fullfile(outdir, ['decoherence_rates_icond', num2str(icond), '.txt']);
    r_ij = importdata(filename);
    rates = complex(r_ij, zeros(nst, nst));
end
for n = 1:num_sh_traj
    me_es{1} = ElectronicStructure_set_state(me_es{1}, init_state);
    me_es{1}.t_m(1) = 0.0;
    for i = 0:(sz - 1)
        if i > 0
            me_es{i + 1} = ElectronicStructure_insertion(me_es{i + 1}, me_es{i});
        end
        me_es{i + 1} = ElectronicStructure_init_hop_prob1(me_es{i + 1});
        [me_es, rates] = propagate_electronic(me_es, i, rates);
        me_es{i + 1} = ElectronicStructure_update_populations(me_es{i + 1});
        if decoherence_method == 0
            [me_es{i + 1}.g, me_es{i + 1}.curr_state] = hop(me_es{i + 1}.g, me_es{i + 1}.curr_state, nst);
        elseif decoherence_method == 1
            me_es{i + 1} = ElectronicStructure_dish_decoherence(me_es{i + 1}, nucl_dt, 1, Temp, rates);
        elseif decoherence_method == 2
            [me_es{i + 1}.g, me_es{i + 1}.curr_state] = hop(me_es{i + 1}.g, me_es{i + 1}.curr_state, nst);
            curr_state = me_es{i + 1}.curr_state;
            me_es{i + 1} = ElectronicStructure_sdm_decoherence(me_es{i + 1}, nucl_dt, curr_state, rates);
            me_es{i + 1} = ElectronicStructure_update_populations(me_es{i + 1});
        elseif decoherence_method == 3
            me_es{i + 1} = ElectronicStructure_gfsh_decohere(me_es{i + 1}, nucl_dt, Temp);
        end
        curr_state = me_es{i + 1}.curr_state;
        sh_pops(i + 1, curr_state + 1) = sh_pops(i + 1, curr_state + 1) + 1;
        se_pops(i + 1, :) = se_pops(i + 1, :) + real(diag(me_es{i + 1}.A)).';
    end
end
se_pops2 = se_pops / num_sh_traj;
sh_pops2 = sh_pops / num_sh_traj;
outfile1 = fullfile(outdir, ['me_pop', num2str(icond)]);
out1 = fopen(outfile1, 'w');
outfile2 = fullfile(outdir, ['out', num2str(icond)]);
out2 = fopen(outfile2, 'w');
for i = 0:(sz - 1)
    fprintf(out1, 'time %d ', i);
    tot = sum(se_pops2(i + 1, :));
    for j = 0:(nst - 1)
        fprintf(out1, 'P(%d)= %.10f ', j, se_pops2(i + 1, j + 1));
    end
    fprintf(out1, 'Total= %.10f\n', tot);
    fprintf(out2, 'time %d ', i);
    for j = 0:(nst - 1)
        fprintf(out2, 'P(%d)= %.10f ', j, sh_pops2(i + 1, j + 1));
    end
    fprintf(out2, '\n');
end
fclose(out1);
fclose(out2);
end
