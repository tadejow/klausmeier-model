#%%
# ==============================================================================
# KOMÓRKA 1: Import pakietów i definicje operatorów
# ==============================================================================
using BifurcationKit, LinearAlgebra, Plots, SparseArrays, Setfield, Statistics
const BK = BifurcationKit

function Laplacian1D_Dirichlet(N, L)
    hx = 2 * L / (N - 1)
    D2 = spdiagm(0 => -2 * ones(N), 1 => ones(N - 1), -1 => ones(N - 1)) / hx^2
    return D2
end

function IntegralOperator(x)
    N = length(x)
    hx = x[2] - x[1]
    weights = ones(N); weights[1] = 0.5; weights[end] = 0.5; weights .*= hx
    kernel(x, y) = (1 / sqrt(pi)) * exp(-(x - y)^2)
    K_matrix = [kernel(xi, xj) for xi in x, xj in x]
    return sparse(K_matrix * diagm(weights))
end

println("Komórka 1 wykonana: Pakiety i operatory zdefiniowane.")

#%%
# ==============================================================================
# KOMÓRKA 2: Definicja modelu matematycznego (Funkcja F i Jakobian J)
# ==============================================================================
function F_system!(f, x, p)
    N = p.N
    u, v = (@view x[1:N]), (@view x[N+1:2N])
    fu, fv = (@view f[1:N]), (@view f[N+1:2N])
    nl_term = u .* u .* v
    mul!(fu, p.Op_u, u); fu .= p.d_u .* fu .- p.B .* u .+ nl_term
    mul!(fv, p.D2, v);   fv .= p.d_v .* fv .- v .+ p.A .- nl_term
    fu[1] = fu[end] = 0.0
    fv[1] = fv[end] = 0.0
    return f
end

function J_system(x, p)
    N = p.N
    u, v = (@view x[1:N]), (@view x[N+1:2N])
    J_uu = p.d_u .* p.Op_u - p.B .* I + spdiagm(0 => 2 .* u .* v)
    J_uv = spdiagm(0 => u .* u)
    J_vu = spdiagm(0 => -2 .* u .* v)
    J_vv = p.d_v .* p.D2 - I - spdiagm(0 => u .* u)
    J = [J_uu J_uv; J_vu J_vv]
    J[1, :] .= 0; J[1, 1] = 1.0
    J[N, :] .= 0; J[N, N] = 1.0
    J[N+1, :] .= 0; J[N+1, N+1] = 1.0
    J[2*N, :] .= 0; J[2*N, 2*N] = 1.0
    return J
end

println("Komórka 2 wykonana: Model matematyczny zdefiniowany.")


#%%
# ==============================================================================
# KOMÓRKA 3: Główna pętla obliczeniowa
# ==============================================================================
B = 0.45; d_u, d_v = 2.0, 80.0
A_min, A_max = 0.1, 1.2
L_values = [1, 5, 10] 
all_diagrams = Dict()

for L in L_values
    println("\n--- Obliczanie dla L = $L ---")
    N = 60
    domain = LinRange(-L, L, N)
    D2_laplace = Laplacian1D_Dirichlet(N, L)
    K_integral = IntegralOperator(domain)
    Op_nonlocal = K_integral - I
    lens = @optic _.A
    params_nonlocal = (N=N, L=L, B=B, d_u=d_u, d_v=d_v, A=A_max, Op_u=Op_nonlocal, D2=D2_laplace)
    params_local = (N=N, L=L, B=B, d_u=d_u, d_v=d_v, A=A_max, Op_u=D2_laplace, D2=D2_laplace)
    record_sol(x, p; kwargs...) = (u_max=maximum(view(x, 1:N)), u_avg=mean(view(x, 1:N)))
    opt_newton = NewtonPar(tol = 1e-5, max_iterations = 100, verbose = false)
    opts_br = ContinuationPar(p_min = A_min, p_max = A_max, ds = -0.01, dsmax = 0.05,
        nev = 10, detect_bifurcation = 3, n_inversion = 6, max_steps=2000, newton_options = opt_newton)
    A_start = A_max
    u_hom = (A_start + sqrt(A_start^2 - 4*B^2)) / (2*B)
    v_hom = B / u_hom
    cos_profile = cos.((π/2) .* domain ./ L)
    u0 = u_hom .* cos_profile .+ 0.01 * (rand(N) .- 0.5) .* cos_profile
    v0 = v_hom .* cos_profile
    x0_guess = vcat(u0, v0)
    
    println("Obliczanie diagramu dla modelu nielokalnego...")
    prob_nl = BifurcationProblem(F_system!, x0_guess, params_nonlocal, lens; J=J_system, record_from_solution=record_sol)
    diagram_nl = @time bifurcationdiagram(prob_nl, PALC(), 2, (args...) -> setproperties(opts_br; ds = -0.01, dsmax = 0.02))
    
    println("Obliczanie diagramu dla modelu lokalnego...")
    prob_l = BifurcationProblem(F_system!, x0_guess, params_local, lens; J=J_system, record_from_solution=record_sol)
    diagram_l = @time bifurcationdiagram(prob_l, PALC(), 2, (args...) -> setproperties(opts_br; ds = -0.01, dsmax = 0.02))
    
    all_diagrams[L] = Dict("nonlocal" => diagram_nl, "local" => diagram_l, "domain" => domain, "N" => N)
end
println("\n--- Wszystkie obliczenia zakończone. Gotowy do rysowania. ---")

#%%
# ==============================================================================
# KOMÓRKA 4: Rysowanie wyników (FINALNA, DZIAŁAJĄCA WERSJA)
# ==============================================================================
plot_array = []
L_values = [1, 5, 10] 

# --- Wiersz 1: Diagramy bifurkacyjne ---
for L in L_values
    p = plot(title="L = $L", xlabel="Rainfall A", xlims=(A_max, A_min))
    if haskey(all_diagrams, L)
        
        print("Printujemy diagramy dla L \n", all_diagrams[L])
        print("Printujemy diagramy dla nl L \n", all_diagrams[L]["nonlocal"])
        print("Printujemy diagramy dla loc L \n", all_diagrams[L]["local"])

        diagram_nl = all_diagrams[L]["nonlocal"]
        diagram_l = all_diagrams[L]["local"]
        
        plot!(p, diagram_nl; vars=(:param, :u_max), legend=:topright, label="Max (Non-local)", color=:yellow)
        plot!(p, diagram_nl; vars=(:param, :u_avg), label="Avg (Non-local)", color=:green)
        plot!(p, diagram_l; vars=(:param, :u_max), label="Max (Local)", color=:cyan, line=:dash)
        plot!(p, diagram_l; vars=(:param, :u_avg), label="Avg (Local)", color=:blue, line=:dash)
    end
    A_curve = LinRange(A_min, A_max, 200)
    plot!(p, A_curve, B ./ A_curve, color=:red, lw=2, label="B/A")
    if L == L_values[1]; ylabel!("Biomass (Avg / Max)"); else; plot!(p, yformatter=_->""); end
    push!(plot_array, p)
end

final_plot = plot(plot_array..., layout=(4, 3), size=(1500, 1800), legend_font_pointsize=8,
                  plot_title="Comparison of Bifurcation Diagrams (Julia / BifurcationKit.jl)")
display(final_plot)