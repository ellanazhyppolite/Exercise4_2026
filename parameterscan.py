#%%

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import glob

# -----------------------------------------------------------------------
# Parameter scan script for the electrostatics exercise.
#
# Default use: convergence study of phi(0) vs N (= N1 = N2) in the
# trivial (uniform) case.  Change 'paramstr' and 'variable_array' to
# scan any other parameter.
# -----------------------------------------------------------------------

# Path to compiled executable (adjust if needed)
repertoire     = '/Users/ella/Documents/EPFL/BA4/PHYSNUM/Exercise4_2026/Exercise4_2026/'
executable     = 'engine'
input_filename = 'trivial.in'   # base configuration file

# Base parameters (values here are overwritten by the scan below)
input_parameters = {
    'b'      : 0.05,   # Inner radius [m]
    'R'      : 0.1,    # Outer radius [m]
    'V0'     : 0,      # Boundary potential [V]
    'a0'     : 1,      # Charge density scale [V/m^2]  (unused when trivial=true)
    'dx'     : 1.0/64.0,
    'trivial': 'true', # true: uniform test case
    'N1'     : 5,      # Intervals in [0, b]
    'N2'     : 5,      # Intervals in [b, R]
}

# -----------------------------------------------------------------------
# Choose the parameter to scan
# -----------------------------------------------------------------------
#paramstr       = 'N1'                        # parameter name in engine
#variable_array = 2**np.arange(1, 9)          # N = 2, 4, 8, ..., 256

paramstr       = 'dx'                        # parameter name in engine
#variable_array = [1.0/64, 1.0/32, 1.0/16]  # rajouter des valeurs pour faire la conv après
variable_array = [1.0/64.0]

# Build a label for output directories / filenames
outstr = (f"electrostatics_b_{input_parameters['b']:.2g}"
          f"_R_{input_parameters['R']:.2g}"
          f"_trivial_{input_parameters['trivial']}")




# -----------------------------------------------------------------------
# Create output directory
# -----------------------------------------------------------------------
outdir = f"Scan_{paramstr}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)

# -----------------------------------------------------------------------
# Run the scan
# -----------------------------------------------------------------------
for val in variable_array:

    params = input_parameters.copy()
    params[paramstr] = val
    # For a convergence study keep N1 = N2
    if paramstr == 'N1':
        params['N2'] = val

    output_file = f"{outstr}_{paramstr}_{val}"
    output_path = os.path.join(outdir, output_file)

    # Build the command-line parameter string
    param_string = " ".join(f"{k}={v}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string} output={output_path}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print("Done.")




#%%






# ============================================================
# Read tir output files  (*_tir.out)
# ============================================================

files = sorted(glob.glob(os.path.join(outdir, "*_tir.out")))

if len(files) == 0:
    print("No output files found! Check that the C++ code ran correctly.")

datasets = []
param_values = []

for f in files:
    name = os.path.basename(f)      
    name = name.replace("_tir.out", "")
    par ts = name.split("_")


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
# Plot: alpha vs dx 
# ============================================================

alphas = []
phis   = []
Es   = []

for data in datasets:

    alphas.append(data[:, 0])
    phis.append(data[:, 1])
    Es.append(data[:, 2])


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
fig.savefig(os.path.join(fig_dir, "tir_alpha_vs_dx.png"), dpi=150)
plt.show()

# ============================================================
# Plot: error vs dx 
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