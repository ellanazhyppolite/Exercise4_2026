#%%

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import subprocess

# Parameters
repertoire = '/Users/ella/Documents/EPFL/BA4/PHYSNUM/Exercise4_2026'
executable = '/engine'
input_filename = 'configuration.in.example'

input_parameters = {
    # Physical parameters
    'b'       : 0.05,
    'R'       : 0.1,
    'V0'      : 0,
    'a0'      : 1,
    'trivial' : 'true',  

    # Numerical parameters
    'N1'      : 5,
    'N2'      : 5,
    'dx'      : 1.0/64,  

    'output'  : 'trivial'
}

# -------------------------------------------------
paramstr = 'dx'
variable_array = [1.0/64, 1.0/32, 1.0/16] 

outstr = f"exercise4_q43a"

# -------------------------------------------------
# Create output directory
# -------------------------------------------------
outdir = outstr
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)

for i in range(len(variable_array)):

    params = input_parameters.copy()
    params[paramstr] = variable_array[i]


    output_prefix = os.path.join(outdir, f"{outstr}_{paramstr}_{variable_array[i]:.6g}")
    params['output'] = output_prefix

    param_string = " ".join(f"{k}={v}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string}"
    )

    print("Running:", cmd)
    subprocess.run(cmd, shell=True)
    print("Done.")


# ============================================================
# Read shooting output files  (*_shooting.out)
# ============================================================

files = sorted(glob.glob(os.path.join(outdir, "*_shooting.out")))

if len(files) == 0:
    print("No output files found! Check that the C++ code ran correctly.")

datasets = []
param_values = []

for f in files:
    name = os.path.basename(f)      
    name = name.replace("_shooting.out", "")
    parts = name.split("_")


    try:
        value = float(parts[-1])
    except ValueError:
        print(f"Could not parse value from filename: {name}")
        continue

    try:
        data = np.loadtxt(f)
    except Exception as e:
        print(f"Could not load {f}: {e}")
        continue

    # make sure data is 2D even if only one row
    if data.ndim == 1:
        data = data[np.newaxis, :]

    datasets.append(data)
    param_values.append(value)

print(f"Found {len(datasets)} datasets.")

if len(datasets) == 0:
    raise RuntimeError("No datasets loaded, stopping here.")

# Sort by parameter value
order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]

# ============================================================
# Figures directory
# ============================================================
fig_dir = os.path.join(outdir, "figures")
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Plot: alpha found vs dx  (convergence of shooting method)
# ============================================================

xs     = []
alphas = []
phis   = []
Es   = []

for data in datasets:

    xs.append(data[:, 0])
    alphas.append(data[:, 1])
    phis.append(data[:, 2])
    Es.append(data[:, 3])


R    = input_parameters['R']
a0   = input_parameters['a0']
exact_alpha = a0 * R / np.pi

fig, ax = plt.subplots()
ax.plot(param_values, alphas, 'o-', label='Numerical $E_r(R)$')
ax.axhline(exact_alpha, color='r', linestyle='--', label=f'Exact = {exact_alpha:.4f}')
ax.set_xlabel('dx')
ax.set_ylabel('$E_r(R)$')
ax.set_title('Shooting method: found $E_r(R)$ vs step size')
ax.legend()
fig.savefig(os.path.join(fig_dir, "shooting_alpha_vs_dx.png"), dpi=150)
plt.show()

# ============================================================
# Plot: error vs dx  (log-log, convergence order)
# ============================================================

errors = np.abs(np.array(alphas) - exact_alpha)

fig, ax = plt.subplots()
ax.loglog(param_values, errors, 'o-', label='Error')
# reference line for order 4 (RK4)
dx_ref = np.array(param_values)
ax.loglog(dx_ref, errors[0] * (dx_ref / dx_ref[0])**4, '--', label='$O(dx^4)$')
ax.set_xlabel('dx')
ax.set_ylabel('$|E_r(R)^{num} - E_r(R)^{exact}|$')
ax.set_title('Convergence of shooting method (RK4)')
ax.legend()
fig.savefig(os.path.join(fig_dir, "shooting_convergence.png"), dpi=150)
plt.show()

# %%