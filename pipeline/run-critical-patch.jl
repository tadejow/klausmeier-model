using LinearAlgebra
using Plots
using ProgressMeter
using LaTeXStrings # For nice plot labels

#=
The IntegralOperator struct and its constructor pre-compute the
discretized kernel and integration weights for a given spatial grid.
=#
struct IntegralOperator
    kernel::Matrix{Float64}
    weights::Vector{Float64}
end

function IntegralOperator(x::AbstractVector)
    dx = x[2] - x[1]
    
    # Equivalent to np.meshgrid(x, x, indexing='ij')
    grid_x = [i for i in x, j in x]
    grid_y = [j for i in x, j in x]
    
    kernel = exp.(-(grid_x .- grid_y).^2)
    
    weights = ones(length(x))
    weights[begin] *= 0.5 # Trapezoidal rule weights
    weights[end] *= 0.5
    weights .*= dx
    
    return IntegralOperator(kernel, weights)
end

"Applies the integral operator to a function u."
function compute(op::IntegralOperator, u::AbstractVector)
    return (1 / sqrt(pi)) * (op.kernel * (u .* op.weights))
end

#=
The LaplacianFD struct and its constructor pre-compute the
finite-difference matrix for the 1D Laplacian operator.
=#
struct LaplacianFD
    D2::Matrix{Float64}
end

function LaplacianFD(x::AbstractVector)
    N = length(x)
    dx = x[2] - x[1]
    
    # Create a tridiagonal matrix using diagm from LinearAlgebra
    D2 = diagm(
        0 => -2 * ones(N),
        1 => ones(N - 1),
       -1 => ones(N - 1)
    ) / dx^2
    
    return LaplacianFD(D2)
end

"Applies the Laplacian operator to a function u."
function compute(op::LaplacianFD, u::AbstractVector)
    return op.D2 * u
end

# --- Simulation Helper Function ---
function run_simulation(L, N, u_init, v_init; model_type="integral")
    domain_omega = range(-L, L, length=N)

    integral_operator = IntegralOperator(domain_omega)
    laplacian_fd = LaplacianFD(domain_omega)

    u_old = copy(u_init)
    v_old = copy(v_init)

    # Enforce zero Dirichlet boundary conditions
    u_old[begin] = u_old[end] = 0
    v_old[begin] = v_old[end] = 0
    
    # Main simulation loop
    for j in 1:500000 # Use a for loop with a max iteration count
        if model_type == "integral"
            u_spatial_term = d_u * compute(integral_operator, u_old) - d_u * u_old
        elseif model_type == "laplacian"
            u_spatial_term = (d_u / 2) * compute(laplacian_fd, u_old)
        else
            error("Unknown model type. Use 'integral' or 'laplacian'.")
        end

        u_new = u_old .+ ht .* (u_spatial_term .+ u_old.^2 .* v_old .- B .* u_old)
        v_new = v_old .+ ht .* (d_v * compute(laplacian_fd, v_old) .- u_old.^2 .* v_old .- v_old .+ A)

        u_new[begin] = u_new[end] = 0
        v_new[begin] = v_new[end] = 0

        if norm(u_new - u_old) < tol
            biomass = sum(u_new) / length(u_new)
            vegetation = hcat(collect(domain_omega), u_new)
            return biomass, vegetation
        end

        if any(isnan, u_new) || any(isnan, v_new)
            return NaN, nothing
        end

        u_old, v_old = copy(u_new), copy(v_new)
    end
    
    # If loop finishes, it failed to converge
    return NaN, nothing
end

# --- Main Simulation ---
println("Starting simulations... (this may take a few minutes)")

# Parameters
const A, B, d_u, d_v = 1.8, 0.45, 5.0, 0.7
const ht = 0.0001
const K = 50
const L_min, L_max = 1.0, 100.0
const tol = 1e-5

# Julia doesn't have geomspace, so we create it manually
L_vals = exp.(range(log(L_min), log(L_max), length=K))

# Pre-allocate result containers with specific types for performance
biomasses_integral = zeros(K)
vegetations_integral = Vector{Union{Matrix{Float64}, Nothing}}(undef, K)
biomasses_laplacian = zeros(K)
vegetations_laplacian = Vector{Union{Matrix{Float64}, Nothing}}(undef, K)

@showprogress "Simulating across L values: " for (i, L) in enumerate(L_vals)
    N = 10
    for j in 1:100
      if L > j 
        N = (j+1) * 10
      end
    end

    u_init = ones(N) .* ((A + sqrt(A^2 - 4*B^2)) / (2 * B))
    v_init = ones(N) .* ((2 * B^2) / (A + sqrt(A^2 - 4*B^2)))

    biomass_int, veg_int = run_simulation(L, N, u_init, v_init, model_type="integral")
    biomasses_integral[i] = biomass_int
    vegetations_integral[i] = veg_int

    biomass_lap, veg_lap = run_simulation(L, N, u_init, v_init, model_type="laplacian")
    biomasses_laplacian[i] = biomass_lap
    vegetations_laplacian[i] = veg_lap
end

# --- Analysis and Plotting ---
uniform_solution_val = (A + sqrt(A^2 - 4*B^2)) / (2*B)
zero_threshold = 0.1
f_size = 12 # Adjusted font size for Plots.jl rendering

# Find critical L values
L_crit_integral, L_crit_laplacian = NaN, NaN

near_zero_indices_int = findall(b -> !isnan(b) && b < zero_threshold, biomasses_integral)
if !isempty(near_zero_indices_int)
    L_crit_integral = L_vals[last(near_zero_indices_int)]
end

near_zero_indices_lap = findall(b -> !isnan(b) && b < zero_threshold, biomasses_laplacian)
if !isempty(near_zero_indices_lap)
    L_crit_laplacian = L_vals[last(near_zero_indices_lap)]
end

println("Critical L (Integral / Non-local):  $(round(L_crit_integral, digits=4))")
println("Critical L (Laplacian / Local): $(round(L_crit_laplacian, digits=4))")

# --- Plotting with Plots.jl ---

# Define a layout: left plot takes 60% width, right side is a 3x3 grid
plot_layout = @layout [a{0.6w} grid(3,3)]

# Create the main plot (plot 'a' in the layout)
p_main = plot(L_vals, biomasses_integral,
    title="Convergence to spatially uniform biomass density",
    xaxis=("Patch half-width L", :log),
    yaxis="Average biomass density",
    label="Non-local model (integral)",
    seriestype=:scatter, color=:green, markersize=3,
    legend=:bottomright, gridalpha=0.4,
    titlefontsize=f_size+2, tickfontsize=f_size-2, legendfontsize=f_size-1, guidefontsize=f_size
)
scatter!(p_main, L_vals, biomasses_laplacian, label="Local model (Laplacian)", color=:blue, markersize=3, markershape=:square)
hline!(p_main, [uniform_solution_val], label=L"Spatially uniform solution $v_{*, 3}$", color=:orange, linestyle=:dash, lw=2)

# Add vertical lines for critical L
if !isnan(L_crit_integral)
    vline!(p_main, [L_crit_integral], label="", color=:green, linestyle=:dash, lw=2)
    annotate!(p_main, L_crit_integral * 1.1, uniform_solution_val*0.9, text("L_crit ≈ $(round(L_crit_integral, digits=2))", :green, :left, f_size-2))
end
if !isnan(L_crit_laplacian)
    vline!(p_main, [L_crit_laplacian], label="", color=:blue, linestyle=:dash, lw=2)
    annotate!(p_main, L_crit_laplacian * 0.9, uniform_solution_val*0.9, text("L_crit ≈ $(round(L_crit_laplacian, digits=2))", :blue, :right, f_size-2))
end

# Create the side plots (the 3x3 grid)
side_plots = []
L_targets = [0.5, 1.5, 2.5, 3.5, 5.0, 10.0, 20.0, 50.0, 100.0]

for (i, target_L) in enumerate(L_targets)
    # Find index of L_vals closest to the target L
    _, idx = findmin(abs.(L_vals .- target_L))
    actual_L = L_vals[idx]
    
    # Create an empty plot for this subplot position
    sp = plot(title="L = $(round(actual_L, digits=2))", 
              ylims=(-0.1, uniform_solution_val * 1.1), 
              gridalpha=0.4,
              legend=false,
              titlefontsize=f_size, tickfontsize=f_size-2
    )

    # Plot local (Laplacian) profile
    if vegetations_laplacian[idx] !== nothing
        data = vegetations_laplacian[idx]
        scatter!(sp, data[:, 1], data[:, 2], color=:blue, markershape=:square, markersize=2, label="Local")
    end

    # Plot non-local (integral) profile
    if vegetations_integral[idx] !== nothing
        data = vegetations_integral[idx]
        scatter!(sp, data[:, 1], data[:, 2], color=:green, markershape=:circle, markersize=2, label="Non-local")
    end
    
    hline!(sp, [uniform_solution_val], color=:orange, linestyle=:dash, lw=2, label="")

    # Add axis labels to specific subplots to avoid clutter
    if i == 4 # Middle left
        ylabel!(sp, "Biomass density", guidefontsize=f_size)
    end
    if i == 8 # Bottom center
        xlabel!(sp, "x", guidefontsize=f_size)
    end

    push!(side_plots, sp)
end

# Combine all plots into the final figure
plot(p_main, side_plots..., layout=plot_layout, size=(1600, 900), margin=5Plots.mm)