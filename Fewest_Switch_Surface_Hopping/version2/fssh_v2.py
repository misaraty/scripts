import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from scipy.linalg import eigh
from findiff import FinDiff

"""
Tully Model 1 Simulation with Surface Hopping Algorithm
======================================================

This script performs the following:
1. Defines two diabatic Hamiltonians (Vplus and Vminus) depending on position x.
2. Computes adiabatic eigenstates and eigenvalues.
3. Calculates the nonadiabatic coupling term d12(x).
4. Saves energy levels and coupling to files.
5. Interpolates these data using splines.
6. Simulates nuclear motion using a surface hopping method.
7. Plots potential surfaces, coupling, and hopping probabilities.
"""

# Set working directory to script location
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Constants from Tully model
A0 = 0.01
B0 = 1.6
C0 = 0.005
D0 = 1.0

# Define diabatic Hamiltonians

def Vplus(x):
    return np.array([
        [A0 * (1 - np.exp(-B0 * x)), C0 * np.exp(-D0 * x**2)],
        [C0 * np.exp(-D0 * x**2), -A0 * (1 - np.exp(-B0 * x))]
    ])

def Vminus(x):
    return np.array([
        [-A0 * (1 - np.exp(B0 * x)), C0 * np.exp(-D0 * x**2)],
        [C0 * np.exp(-D0 * x**2), A0 * (1 - np.exp(B0 * x))]
    ])

# Compute eigenvalues and eigenvectors

def compute_eigen(x):
    H = Vplus(x) if x > 0 else Vminus(x)
    eps, vecs = eigh(H)
    vecs = vecs / np.linalg.norm(vecs, axis=0)
    return eps, vecs

def epsilon1(x):
    return compute_eigen(x)[0][0]

def epsilon2(x):
    return compute_eigen(x)[0][1]

# Compute nonadiabatic coupling d12(x)

def d12(x, dx=1e-3):
    x_grid = np.linspace(x - 2*dx, x + 2*dx, 5)
    phi2_vals = np.array([compute_eigen(x_)[1][:, 1] for x_ in x_grid])
    diff_op = FinDiff(0, dx, 1)
    dphi2_dx = diff_op(phi2_vals)[2]  # central derivative
    phi1 = compute_eigen(x)[1][:, 0]
    return np.dot(phi1.conj(), dphi2_dx).real

# Save ε1(x), ε2(x), d12(x)

xgrid = np.linspace(-10, 10, 100)
np.savetxt("Tully-model1-e1.dat", np.array([[x, epsilon1(x)] for x in xgrid]))
np.savetxt("Tully-model1-e2.dat", np.array([[x, epsilon2(x)] for x in xgrid]))
np.savetxt("Tully-model1-d12.dat", np.array([[x, d12(x)] for x in xgrid]))

# Load data and spline fit
xgrid_data, e1_data = np.loadtxt("Tully-model1-e1.dat").T
_, e2_data = np.loadtxt("Tully-model1-e2.dat").T
_, d12_data = np.loadtxt("Tully-model1-d12.dat").T

e1fit = interp.splrep(xgrid_data, e1_data)
e2fit = interp.splrep(xgrid_data, e2_data)
d12fit = interp.splrep(xgrid_data, d12_data)

e1 = lambda x: interp.splev(x, e1fit)
e2 = lambda x: interp.splev(x, e2fit)
d12 = lambda x: interp.splev(x, d12fit)

# Plotting ε1(x), ε2(x), and d12(x)
xfine = np.linspace(-10, 10, 1000)
fig, ax = plt.subplots()
ax.set_xlabel("Ion Position (x)", fontsize=14)
ax.set_ylabel("Potential and Coupling", fontsize=14)
ax.plot(xfine, e1(xfine), "-", lw=2, label="$\\epsilon_1$")
ax.plot(xfine, e2(xfine), "--", lw=2, label="$\\epsilon_2$")
ax.plot(xfine, -d12(xfine)/50., "-.", lw=2, label="$-d_{12}/50$")
ax.legend(fontsize=14, loc="upper right")
fig.tight_layout()
plt.savefig("single.jpg", dpi=300)

# Derivatives for forces
e1der = lambda x: interp.splev(x, interp.splder(e1fit))
e2der = lambda x: interp.splev(x, interp.splder(e2fit))

# Surface hopping trajectory simulation

pos_min, pos_max = -5.0, 5.0

def run1(x0, v0, c1, nstep, dt=0.1, m=2000.):
    hbar = 1.0
    c2 = np.sqrt(1. - c1*c1)
    a11, a12, a21, a22 = c1*np.conj(c1), c1*np.conj(c2), np.conj(c1*np.conj(c2)), c2*np.conj(c2)
    pos, vel = x0, v0

    pos_trace, vel_trace = np.zeros(nstep), np.zeros(nstep)
    a11_trace, a22_trace = np.zeros(nstep), np.zeros(nstep)
    energy_trace, pot_trace, prob_trace = np.zeros(nstep), np.zeros(nstep), np.zeros(nstep)

    on_surface1 = True
    surface, surface_der = e1, e1der

    for t in range(nstep):
        pos_trace[t], vel_trace[t] = pos, vel
        pot_trace[t] = surface(pos)
        energy_trace[t] = 0.5 * m * vel**2 + surface(pos)
        a11_trace[t], a22_trace[t], prob_trace[t] = a11.real, a22.real, 0.0

        # Propagate coefficients
        c1 += dt * (c1 * e1(pos) / (1j * hbar) - c2 * vel * d12(pos))
        c2 += dt * (c2 * e2(pos) / (1j * hbar) + c1 * vel * d12(pos))

        a11, a12, a21, a22 = c1*np.conj(c1), c1*np.conj(c2), np.conj(c1*np.conj(c2)), c2*np.conj(c2)
        b12 = -2. * (np.conj(a12) * vel * d12(pos)).real
        b21 =  2. * (np.conj(a21) * vel * d12(pos)).real
        switch_prob = dt * (b12 / a11 if on_surface1 else b21 / a11)
        prob_trace[t] = switch_prob.real

        if switch_prob.real > np.random.rand():
            new_kinetic = energy_trace[t] - (e2(pos) if on_surface1 else e1(pos))
            if new_kinetic >= 0:
                on_surface1 = not on_surface1
                surface = e2 if not on_surface1 else e1
                surface_der = e2der if not on_surface1 else e1der
                vel = np.sign(vel) * np.sqrt(2. * new_kinetic / m)

        a0 = -surface_der(pos) / m
        pos += vel * dt + 0.5 * a0 * dt**2
        a1 = -surface_der(pos) / m
        vel += 0.5 * (a0 + a1) * dt

        if pos >= pos_max or pos <= pos_min:
            break

    transmission = pos >= pos_max
    reflection = pos <= pos_min

    return pos_trace[:t], vel_trace[:t], pot_trace[:t], energy_trace[:t], a11_trace[:t], a22_trace[:t], t, on_surface1, transmission, reflection, prob_trace[:t]

# Run multiple trajectories and collect statistics
klist = [2.5, 4.5, 10, 15, 20, 25, 30, 35]
ntraj = 200
ntrans1 = np.zeros(len(klist))
nref1   = np.zeros(len(klist))
ntrans2 = np.zeros(len(klist))
nref2   = np.zeros(len(klist))
not_finished = np.zeros(len(klist))

for i, k in enumerate(klist):
    for _ in range(ntraj):
        vo = k / 2000.
        result = run1(-5.0, vo, 1.0, 4000, dt=2.0)
        _, _, _, _, _, _, _, on1, trans, ref, _ = result

        if on1 and trans:
            ntrans1[i] += 1
        elif on1 and ref:
            nref1[i] += 1
        elif not on1 and trans:
            ntrans2[i] += 1
        elif not on1 and ref:
            nref2[i] += 1
        else:
            not_finished[i] += 1

# Plot statistics
fig = plt.figure()
plt.subplots_adjust(hspace=0.001)

ax1 = fig.add_subplot(311)
ax1.plot(klist, ntrans1 / ntraj, "o-")
ax1.set_ylabel("trans. on 1", fontsize=14)
ax1.set_ylim(-0.2, 1.2)
ax1.set_yticks([0, 0.5, 1])

ax2 = fig.add_subplot(312, sharex=ax1)
ax2.plot(klist, nref1 / ntraj, "o-")
ax2.set_ylabel("ref. on 1", fontsize=14)
ax2.set_ylim(-0.2, 1.2)
ax2.set_yticks([0, 0.5, 1])

ax3 = fig.add_subplot(313, sharex=ax1)
ax3.plot(klist, ntrans2 / ntraj, "o-")
ax3.set_ylabel("trans. on 2", fontsize=14)
ax3.set_xlabel("$k_o$ (a.u.)", fontsize=14)
ax3.set_ylim(-0.2, 1.2)
ax3.set_yticks([0, 0.5, 1])

plt.savefig("model-1-prob.jpg", dpi=300)