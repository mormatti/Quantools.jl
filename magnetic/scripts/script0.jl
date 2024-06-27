include("../../main.jl")

wd = "/Users/Mattia/Code/QuantumProjects/magnetic"

using LinearAlgebra
using Plots
using .QuantumProjects
using .Utils
using .BasicTools
using .LinearAlgebraExtensions
using .Magnetic

let # Beginning of script

St = [true, false]
Lk = [-1,0,1]
n = 3
Lx = 2
Ly = 2

BasisPhys = []
BasisVort = []

# We do the same loop of before but considering a 3x2 lattice
for s1 in St, s2 in St, s3 in St, s4 in St, l1 in Lk, l2 in Lk, l3 in Lk, l4 in Lk
    C = Config(Ly, Lx, n)
    C.sites[1,1] = s1
    C.sites[1,2] = s2
    C.sites[2,1] = s3
    C.sites[2,2] = s4
    C.xlinks[1,1] = l1
    C.xlinks[1,2] = l2
    C.ylinks[1,1] = l3
    C.ylinks[2,1] = l4
    if isGaugeInvariant(C)
        push!(BasisPhys, deepcopy(C))
    end
    if isGaugeInvariant(C) && withoutCharge(C)
        push!(BasisVort, deepcopy(C))
    end
end

BasisExt = tensorProductBasis(BasisPhys, BasisVort)
Np = length(BasisPhys)
Nv = length(BasisVort)
Ne = length(BasisExt)

println("Number of physical states: ", Np)
println("Number of vortices states: ", Nv)
println("Number of extended states: ", Ne)

global M = matrix(BasisExt, C -> mapping(C))

global P = (M' * M) / Nv
global Π = (M * M') / Nv

Plotting.matrix(M, filepath = "$wd/images/M.png", dpi = 500)
Plotting.matrix(P, filepath = "$wd/images/P.png", dpi = 500)
Plotting.matrix(Π, filepath = "$wd/images/Π.png", dpi = 500)

println()
println("Rank")
println("Rk(M) = ", rank(M))
println("Rk(P) = ", rank(P))
println("Rk(Π) = ", rank(Π))

println()
println("Trace")
println("Tr(M) = ", tr(M))
println("Tr(P) = ", tr(P))
println("Tr(Π) = ", tr(Π))

println()
println("Determinant")
println("det(M)   = ", det(M))
println("det(P) = ", det(P))
println("det(Π) = ", det(Π))

println()
println("Hermiticity")
println("M Hermitian: ", M ≈ M')
println("P Hermitian: ", P ≈ P')
println("Π Hermitian: ", Π ≈ Π')

println()
println("Idempotency")
println("M idempotent: ", M ≈ M^2)
println("P idempotent: ", P ≈ P^2)
println("Π idempotent: ", Π ≈ Π^2)

global Hp = HamiltonianMatrix(BasisPhys, 1, 1, 1, 1)
global He = HamiltonianMatrix(BasisExt, 1, 1, 1, 1)

Plotting.matrix(Hp, filepath = "$wd/images/Hp.png", dpi = 500)
Plotting.matrix(He, filepath = "$wd/images/He.png", dpi = 500)

HpSpec = real(eigvals(Hp))
HeSpec = real(eigvals(He))
ProdSpec = real(eigvals(P * He))

function drawSpec(eigenvals; filepath::String = "spectrum.png", dpi::Integer = 500)
    plot(real(eigenvals), dpi = dpi)
    savefig(filepath)
end

drawSpec(HpSpec, filepath = "$wd/images/Spec-Hp.png")
drawSpec(HeSpec, filepath = "$wd/images/Spec-He.png")
drawSpec(ProdSpec, filepath = "$wd/images/Spec-Prod.png")

indexToEliminate = []
for i in eachindex(ProdSpec)
    if abs(ProdSpec[i]) < 1e-7
        push!(indexToEliminate, i)
    end
end
ProdSpec = deleteat!(ProdSpec, indexToEliminate)

println("Number of non-zero eigenvalues: ", length(ProdSpec))
drawSpec(ProdSpec, filepath = "$wd/images/Spec-ProdCut.png")

println("Computing eigenvectors and eigenvalues...")
HeEigen = eigen(He)
ProdEigen = eigen(P * He)

println("Groundstate energy of He: ", minimum(HeSpec))
println("Groundstate energy of P * He: ", minimum(ProdSpec))
println("Norm of difference of groundstates:", norm(HeEigen.vectors[:,1] - ProdEigen.vectors[:,1]))
println("Norm of difference of second excited state:", norm(HeEigen.vectors[:,2] - ProdEigen.vectors[:,2]))

# end of script
end