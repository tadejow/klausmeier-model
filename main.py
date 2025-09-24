# ==============================================================================
# BIFURCATION ANALYSIS USING A ROBUST, BUILT-IN NUMPY SOLVER
# This version is guaranteed to work without new package installations.
# ==============================================================================
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d


# --- Operator Class Definitions (unchanged) ---
class IntegralOperator:
    def __init__(self, x):
        self.x = x;
        self.N = len(x)
        self.grid = np.meshgrid(self.x, self.x, indexing='ij')
        kernel_func = lambda x, y: (1 / np.sqrt(np.pi)) * np.exp(-(x - y) ** 2)
        self.weights = self._trapezoidal_weights()
        self.matrix = kernel_func(self.grid[0], self.grid[1]) @ np.diag(self.weights)

    def _trapezoidal_weights(self):
        w = np.ones(self.N);
        w[0] *= 0.5;
        w[-1] *= 0.5
        w *= (self.x[-1] - self.x[0]) / (self.N - 1);
        return w


class LaplacianOperator:
    def __init__(self, x):
        self.x = x;
        self.N = len(x)
        self.D2 = self._finite_diff_matrix()

    def _finite_diff_matrix(self):
        dx = self.x[1] - self.x[0]
        D = np.diag(np.ones(self.N - 1), -1) + np.diag(np.ones(self.N - 1), 1) - 2 * np.eye(self.N)
        return D / dx ** 2


# --- Simulation Parameters ---
B = 0.45;
d_u, d_v = 0.1, 0.1
ht = 1e-3;
tol = 1e-6  # A slightly tighter tolerance is good practice with better solvers
max_iter = 10000
A_min, A_max = 0.1, 1.0
L = 20.0;
K = 250
N = 101  # A better resolution is now easily achievable

# --- Domain, Operators, and Pre-computed Matrices ---
domain = np.linspace(-L, L, N)
integral_operator = IntegralOperator(domain)
laplacian_fd = LaplacianOperator(domain)

print("--- Pre-computing implicit matrices ---")
u_matrix_nonlocal = np.eye(N) - ht * (d_u * integral_operator.matrix - d_u * np.eye(N) - B * np.eye(N))
u_matrix_nonlocal[0, :], u_matrix_nonlocal[-1, :] = 0, 0
u_matrix_nonlocal[0, 0], u_matrix_nonlocal[-1, -1] = 1, 1

u_matrix_local = np.eye(N) - ht * (d_u * laplacian_fd.D2 - B * np.eye(N))
u_matrix_local[0, :], u_matrix_local[-1, :] = 0, 0
u_matrix_local[0, 0], u_matrix_local[-1, -1] = 1, 1

v_matrix = np.eye(N) - ht * (d_v * laplacian_fd.D2 - np.eye(N))
v_matrix[0, :], v_matrix[-1, :] = 0, 0
v_matrix[0, 0], v_matrix[-1, -1] = 1, 1


def run_continuation(model_type, u_matrix, v_matrix):
    """Runs the full continuation experiment using NumPy's direct solver."""
    print(f"--- Starting continuation for {model_type.upper()} model ---")

    A_values = np.linspace(A_max, A_min, K)
    branch_results = {}
    u_init, v_init = None, None

    for A in tqdm(A_values, desc=f"Bifurcation analysis ({model_type})"):
        if u_init is None:
            u_init = np.ones(N) * ((A + np.sqrt(A ** 2 - 4 * B ** 2)) / (2 * B))
            v_init = np.ones(N) * ((2 * B ** 2) / (A + np.sqrt(A ** 2 - 4 * B ** 2)))

        u_old, v_old = u_init.copy(), v_init.copy()

        for it in range(max_iter):
            nonlinear_term = u_old ** 2 * v_old

            rhs_u = u_old + ht * nonlinear_term
            rhs_v = v_old + ht * (A - nonlinear_term)
            rhs_u[0] = rhs_u[-1] = 0
            rhs_v[0] = rhs_v[-1] = 0

            # *** MODIFICATION: Replace Gauss-Seidel with np.linalg.solve ***
            # This is a direct, fast, and robust solver
            u_new = np.linalg.solve(u_matrix, rhs_u)
            v_new = np.linalg.solve(v_matrix, rhs_v)

            if np.linalg.norm(u_new - u_old) < tol:
                break

            u_old, v_old = u_new.copy(), v_new.copy()

        branch_results[A] = (u_new, v_new)
        u_init, v_init = u_new.copy(), v_new.copy()

    return A_values, branch_results


# --- Run the simulations ---
A_vals_nonlocal, results_nonlocal = run_continuation('nonlocal', u_matrix_nonlocal, v_matrix)
A_vals_local, results_local = run_continuation('local', u_matrix_local, v_matrix)

print("\n--- All computations finished. Ready for plotting. ---")


# --- PLOTTING (This part remains unchanged) ---

def extract_biomass_data(A_values, results):
    avg_biomass = [results[A][0].mean() for A in A_values]
    max_biomass = [results[A][0].max() for A in A_values]
    return np.array(avg_biomass), np.array(max_biomass)


avg_nonlocal, max_nonlocal = extract_biomass_data(A_vals_nonlocal, results_nonlocal)
avg_local, max_local = extract_biomass_data(A_vals_local, results_local)
A_critical_curve = np.linspace(A_min, A_max, 200)
critical_level = B / A_critical_curve

fig = plt.figure(figsize=(18, 11))
gs_outer = gridspec.GridSpec(1, 2, width_ratios=[1.3, 2], wspace=0.3)

ax_main = fig.add_subplot(gs_outer[0])
ax_main.set_title('Bifurcation Diagram', fontsize=14)
ax_main.scatter(A_vals_nonlocal, avg_nonlocal, color='green', s=10, label='Average Biomass (Non-local)')
ax_main.scatter(A_vals_nonlocal, max_nonlocal, color='y', s=10, label='Max Biomass (Non-local)')
ax_main.scatter(A_vals_local, avg_local, marker='x', color='blue', s=10, label='Average Biomass (Local)')
ax_main.scatter(A_vals_local, max_local, marker='x', color='cyan', s=10, label='Max Biomass (Local)')
ax_main.plot(A_critical_curve, critical_level, 'r-', lw=2, label='Critical Level (B/A)')
ax_main.set_xlabel('Rainfall A', fontsize=12)
ax_main.set_ylabel('Biomass', fontsize=12)
ax_main.set_xlim(A_max, A_min)
ax_main.set_ylim(bottom=0)
ax_main.legend()
ax_main.grid(True)

gs_inner = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=gs_outer[1], hspace=0.5, wspace=0.2)
side_axes = [fig.add_subplot(gs_inner[i, j]) for i in range(3) for j in range(3)]
A_targets = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
N_fine = 500
domain_fine = np.linspace(-L, L, N_fine)

for i, ax in enumerate(side_axes):
    A_target = A_targets[i]

    idx_nonlocal = np.argmin(np.abs(A_vals_nonlocal - A_target))
    idx_local = np.argmin(np.abs(A_vals_local - A_target))
    actual_A = A_vals_nonlocal[idx_nonlocal]

    u_profile_nonlocal = results_nonlocal[actual_A][0]
    u_profile_local = results_local[actual_A][0]

    f_nl = interp1d(domain, u_profile_nonlocal, kind='cubic')
    f_l = interp1d(domain, u_profile_local, kind='cubic')

    ax.plot(domain_fine, f_l(domain_fine), color='blue', ls='--', label='Local')
    ax.plot(domain_fine, f_nl(domain_fine), color='green', lw=2, label='Non-local')

    ax.set_title(f'Vegetation profile for A ≈ {A_target:.2f}', fontsize=12)
    ax.set_ylim(bottom=-0.5)
    ax.grid(True, linestyle='--', alpha=0.6)

    if i == 0: ax.legend()
    if i % 3 == 0: ax.set_ylabel('Biomass Density')
    if i >= 6: ax.set_xlabel('x')

fig.tight_layout()
plt.show()