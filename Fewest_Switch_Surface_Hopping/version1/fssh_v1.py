import os
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

# Set working directory
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Load data from files
xgrid_data, e1_data = np.loadtxt("Tully-model1-e1.dat").T
_, e2_data = np.loadtxt("Tully-model1-e2.dat").T
_, d12_data = np.loadtxt("Tully-model1-d12.dat").T

# Interpolate using spline fits
e1fit = interp.splrep(xgrid_data, e1_data)
e2fit = interp.splrep(xgrid_data, e2_data)
d12fit = interp.splrep(xgrid_data, d12_data)

e1 = lambda x: interp.splev(x, e1fit)
e2 = lambda x: interp.splev(x, e2fit)
d12 = lambda x: interp.splev(x, d12fit)

# Define x range for plotting
pos_min = -10.0
pos_max = 10.0
xfine = np.linspace(pos_min, pos_max, 1000)

# Plot potential energy surfaces and coupling
fig, ax = plt.subplots()
ax.set_xlabel("Ion Position (x)", fontsize=14)
ax.set_ylabel("Potential and Coupling", fontsize=14)
ax.plot(xfine, e1(xfine), "-", lw=2, label="$\\epsilon_1$")
ax.plot(xfine, e2(xfine), "--", lw=2, label="$\\epsilon_2$")
ax.plot(xfine, -d12(xfine) / 50., "-.", lw=2, label="$-d_{12}/50$")
ax.legend(fontsize=14, loc="upper right")
fig.tight_layout()
plt.savefig("single.jpg", dpi=300)

# Compute derivatives of interpolated surfaces
e1der = lambda x: interp.splev(x, interp.splder(e1fit))
e2der = lambda x: interp.splev(x, interp.splder(e2fit))

# Reduce simulation region to avoid spline boundary artifacts
pos_min /= 2.
pos_max /= 2.

def run1(x0, v0, c1, nstep, dt=0.1, m=2000.):
    hbar = 1.0
    c2 = np.sqrt(1. - c1 * c1)
    a11, a12, a21, a22 = c1 * np.conj(c1), c1 * np.conj(c2), np.conj(c1 * np.conj(c2)), c2 * np.conj(c2)
    pos, vel = x0, v0

    pos_trace = np.zeros(nstep)
    vel_trace = np.zeros(nstep)
    a11_trace = np.zeros(nstep)
    a22_trace = np.zeros(nstep)
    energy_trace = np.zeros(nstep)
    pot_trace = np.zeros(nstep)
    prob_trace = np.zeros(nstep)

    on_surface1 = True
    surface = e1
    surface_der = e1der
    switch_prob = 0.0

    for t in range(nstep):
        pos_trace[t] = pos
        vel_trace[t] = vel
        pot_trace[t] = surface(pos)
        energy_trace[t] = 0.5 * m * vel * vel + surface(pos)
        a11_trace[t] = a11.real
        a22_trace[t] = a22.real
        prob_trace[t] = switch_prob.real

        c1 += dt * (c1 * e1(pos) / (1j * hbar) - c2 * vel * d12(pos))
        c2 += dt * (c2 * e2(pos) / (1j * hbar) + c1 * vel * d12(pos))

        a11 = c1 * np.conj(c1)
        a12 = c1 * np.conj(c2)
        a21 = np.conj(a12)
        a22 = c2 * np.conj(c2)

        b12 = -2. * (np.conj(a12) * vel * d12(pos)).real
        b21 =  2. * (np.conj(a21) * vel * d12(pos)).real

        switch_prob = dt * (b12 / a11 if on_surface1 else b21 / a11)

        if switch_prob.real > np.random.rand():
            new_kinetic = energy_trace[t] - (e2(pos) if on_surface1 else e1(pos))
            if new_kinetic >= 0:
                on_surface1 = not on_surface1
                surface = e2 if not on_surface1 else e1
                surface_der = e2der if not on_surface1 else e1der
                vel = np.sign(vel) * np.sqrt(2. * new_kinetic / m)

        a0 = -surface_der(pos) / m
        pos += vel * dt + 0.5 * a0 * dt * dt
        a1 = -surface_der(pos) / m
        vel += 0.5 * (a0 + a1) * dt

        if pos >= pos_max or pos <= pos_min:
            break

    transmission = pos >= pos_max
    reflection = pos <= pos_min

    return pos_trace[:t], vel_trace[:t], pot_trace[:t], energy_trace[:t], a11_trace[:t], a22_trace[:t], t, on_surface1, transmission, reflection, prob_trace[:t]

# Simulate multiple initial conditions and trajectories
klist = [2.5, 4.5, 10, 15, 20, 25, 30, 35]
ntraj = 200  # Reduced from 200 for faster runtime
nstep = 4000  # Reduced steps
ntrans1 = np.zeros(len(klist))
ntrans2 = np.zeros(len(klist))
nref1 = np.zeros(len(klist))
nref2 = np.zeros(len(klist))
not_finished = np.zeros(len(klist))

for i, k in enumerate(klist):
    for _ in range(ntraj):
        vo = k / 2000.
        result = run1(-5.0, vo, 1.0, nstep, dt=2.5)  # dt increased from 2.0 to 2.5
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

# Plot the hopping probabilities
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
