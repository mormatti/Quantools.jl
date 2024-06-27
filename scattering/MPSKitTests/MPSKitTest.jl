using MPSKit, MPSKitModels, TensorKit
using ProgressMeter, Plots # for demonstration purposes

L = 16 # length of the chain
D = 4 # bonddimension
init_state = FiniteMPS(L, ℂ^2, ℂ^D)

g_values = 0:0.1:2
Z = @mpoham sum(σᶻ(){i} for i in vertices(FiniteChain(L)))

M = @showprogress map(g_values) do g
    H = periodic_boundary_conditions(transverse_field_ising(; g=g), L)
    groundstate, environment, δ = find_groundstate(init_state, H; verbose=false)
    return abs(sum(expectation_value(groundstate, Z))) / L
end

scatter(g_values, M, xlabel="g", ylabel="M", label="D=$D", title="Magnetization")

propagator(groundstate, H, δ; verbose=false)