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
#repertoire     = '/Users/ella/Documents/EPFL/BA4/PHYSNUM/Exercise4_2026/Exercise4_2026/'
repertoire     = '/Users/eldidi/Desktop/Exercise4_2026/'
executable     = 'engine'
input_filename = 'nontrivial.in'   # base configuration file

# Base parameters (values here are overwritten by the scan below)
input_parameters = {
    'b'      : 0.02,   # Inner radius [m]
    'R'      : 0.1,    # Outer radius [m]
    'V0'     : 0,      # Boundary potential [V]
    'a0'     : 10000,      # Charge density scale [V/m^2]  (unused when trivial=true)
    'dx'     : 1.0/64.0,
    'trivial': 'false', # true: uniform test case
    'N1'     : 5,      # Intervals in [0, b]
    'N2'     : 160,      # Intervals in [b, R]
}

# -----------------------------------------------------------------------
# Choose the parameter to scan
# -----------------------------------------------------------------------
#paramstr       = 'N1'                        # parameter name in engine
#variable_array = 2**np.arange(1, 9)          # N = 2, 4, 8, ..., 256

paramstr       = 'N1'                        # parameter name in engine
variable_array = [160]
#variable_array = [1.0/64, 1.0/32, 1.0/16]  # rajouter des valeurs pour faire la conv après
#variable_array = [1.0/64.0]

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






files = sorted(glob.glob(os.path.join(outdir, "*_divD_rho.out")))

if len(files) == 0:
    print("No output files found! Check that the C++ code ran correctly.")

datasets = []
param_values = []

for f in files:
    name = os.path.basename(f)      
    name = name.replace("_divD_rho.out", "")
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
# Plot: D vs r + erruer relative 
# ======================= =====================================

r_coords = data[:, 0]        
div_D = data[:, 1]          
rho_sur_eps0 = data[:, 2]


plt.figure(figsize=(10, 6))

plt.plot(r_coords, div_D, 'o', label=r'$\nabla \cdot \mathbf{D} / \epsilon_0$ ', markersize=4)


plt.plot(r_coords, rho_sur_eps0, '-', label=r'$\rho_{lib} / \epsilon_0$ ', alpha=0.7)

plt.xlabel(' r [m]', fontsize=14)
plt.ylabel('divergence et densité de charge [V/m$^2$]', fontsize=14)
plt.legend()

plt.grid(True)


erreur_abs = np.abs(div_D - rho_sur_eps0)  


plt.figure(figsize=(10, 5))


plt.semilogy(r_coords, erreur_abs, 'r-', linewidth=1.5, label=r'Erreur absolue $| \nabla \cdot \mathbf{D}/\epsilon_0 - \rho_{lib}/\epsilon_0 |$')


plt.xlabel('r [m]', fontsize=12)
plt.ylabel('Erreur absolue [$V/m^2$]', fontsize=12)
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.legend()

# Optionnel : marquer la frontière b = 0.02
plt.axvline(x=0.02, color='black', linestyle='--', alpha=0.3, label='Interface b')

plt.tight_layout()
plt.show()




plt.show()