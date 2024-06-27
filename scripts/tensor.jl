include("../utils/main.jl")
using LinearAlgebra
using .Utils.BasicTools
using .Utils.LinearAlgebraExtensions
using .Utils.TensorNetworks
using ITensors
using Plots
using Optim

function opSumIsing(
        L::Int64,
        λ::Float64;
        n::Int64 = 0,
        prefactor::Float64 = 1.0,
        range::Union{Nothing,Int64} = nothing
        )

    @assert L % 2 == 1 # We assert that the number of sites is odd
    jc = (L+1)/2 # The central site position
    χ(j) = j - jc # The distance function
    f(j) = (χ(j))^n
    J, h = 1-λ, λ

    os = OpSum()
    for j ∈ 1:L
        p = prefactor * f(j)
        p = (isnothing(range) ? p : (abs(χ(j)) ≤ range ? p : 0.0))
        if j == 1
            os += p * 0.5 * J, "Sz", 1, "Sz", 2
            os += p * h, "Sx", 1
        elseif j == L
            os += p * 0.5 * J, "Sz", L-1, "Sz", L
            os += p * h, "Sx", L
        else
            os += p * 0.5 * J, "Sz", j-1, "Sz", j
            os += p * 0.5 * J, "Sz", j, "Sz", j+1
            os += p * h, "Sx", j
        end
    end
    
    return os
end

function opSumIsingPBC(
    L::Int64,
    λ::Float64;
    n::Int64 = 0,
    prefactor::Float64 = 1.0,
    range::Union{Nothing,Int64} = nothing
    )

    @assert L % 2 == 1 # We assert that the number of sites is odd
    jc = (L+1)/2 # The central site position
    χ(j) = j - jc # The distance function
    f(j) = (χ(j))^n
    J, h = 1-λ, λ

    os = OpSum()
    for j ∈ 1:L
        p = prefactor * f(j)
        p = (isnothing(range) ? p : (abs(χ(j)) ≤ range ? p : 0.0))
        os += p * J, "Sz", j↻L, "Sz", (j+1)↻L
        os += p * h, "Sx", j↻L
    end

return os
end

function opSumLocalIsing(
        L::Int64, 
        λ::Float64, 
        j::Int64
        )

    @assert L % 2 == 1
    @assert 1 ≤ j ≤ L
    J = 1-λ
    h = λ
    os = OpSum()
    if j == 1
        os += 0.5 * J, "Sz", 1, "Sz", 2
        os += h, "Sx", 1
    elseif j == L
        os += 0.5 * J, "Sz", L-1, "Sz", L
        os += h, "Sx", L
    else
        os += 0.5 * J, "Sz", j-1, "Sz", j
        os += 0.5 * J, "Sz", j, "Sz", j+1
        os += h, "Sx", j
    end
end

function opSumIsingAntisymm(
    L::Int64,
    λ::Float64;
    prefactor::Float64 = 1.0
    )

    @assert L % 2 == 1 # We assert that the number of sites is odd
    jc = (L+1)/2 # The central site position
    χ(j) = j - jc # The distance function
    f(j) = j==0 ? 0 : j>0 ? 1 : -1
    J, h = 1-λ, λ

    os = OpSum()
    for j ∈ 1:L
        p = prefactor * f(j)
        if j == 1
            os += p * 0.5 * J, "Sz", 1, "Sz", 2
            os += p * h, "Sx", 1
        elseif j == L
            os += p * 0.5 * J, "Sz", L-1, "Sz", L
            os += p * h, "Sx", L
        else
            os += p * 0.5 * J, "Sz", j-1, "Sz", j
            os += p * 0.5 * J, "Sz", j, "Sz", j+1
            os += p * h, "Sx", j
        end
    end

    return os
end

function matrixIsing(
    ; L::Int64=11, 
    λ::Float64=0.5
    )

    Sx = [0 0.5; 0.5 0]
    Sz = [0.5 0; 0 -0.5]
    return λ * Σⱼ(Sz, Sz, L) + (1-λ) * Σⱼ(Sx, L)
end

function reflect(ψ::MPS)
    N, Sd = length(ψ), siteinds(ψ)
    ϕ = MPS(N)
    for j in 1:N
        h = N - j + 1
        ϕ[h] = ψ[j] * delta(Sd[j], Sd[h]')
    end
    noprime!(siteinds, ϕ)
    return ϕ
end

function translate(ψ::MPS; dir = "right")
    L = length(ψ)
    ϕ = copy(ψ)
    if dir == "right"
        for j in 1:L-1
            ϕ = swapbondsites(ϕ, j)
        end
    elseif dir == "left"
        for j in L-1:-1:1
            ϕ = swapbondsites(ϕ, j)
        end
    else
        error("Invalid direction")
    end
    return ϕ
end

let
    w = 4
    λ = 0.8
    L = 4 * w + 1
    jc = round(Int64, (L+1)/2)
    nₑ = L # Number of excited levels to be computed, L if all

    # We create the sites and the Hamiltonian
    sites = siteinds("S=1/2", L)
    𝐇 = MPO(opSumIsingPBC(L, λ; prefactor=-1.0), sites)
    𝐡 = [MPO(opSumLocalIsing(L, λ, j), sites) for j ∈ 1:L]

    # DMRG parameters
    χ₀ = 200 # The initial bond dimension of DMRG
    χₘ = 200 # The maximum bond dimension of DMRG
    sₘ = 10000 # The maximum number of sweeps of DMRG

    # We run the DMRG algorithm in order to find the ground state and the first excited states
    observer = DMRGObserver(energy_tol = 10^(-10))
    E₀, 𝛙₀ = dmrg(𝐇, randomMPS(sites; linkdims = χ₀); nsweeps = sₘ, maxdim = χₘ, cutoff = 10^(-15), observer = observer)
    energy_density(𝛟,j) = real(inner(𝛟', 𝐡[j], 𝛟) - inner(𝛙₀', 𝐡[j], 𝛙₀))
    # plot([energy_density(𝛙₀,j) for j ∈ 1:L], label = "|gs⟩")
    parity(𝛟) = inner(𝛟, reflect(𝛟))
    cosk(𝛟) = inner(𝛟', translate(𝛟))
    𝛙::Vector{MPS} = []
    E = []
    for α ∈ 1:nₑ
        Eₙ, 𝛙ₙ = dmrg(𝐇, [𝛙₀; 𝛙], randomMPS(sites; linkdims = χ₀); nsweeps = sₘ, maxdim = χₘ, cutoff = 10^(-15), observer = observer)
        push!(E, Eₙ)
        push!(𝛙, 𝛙ₙ)
        # plot!([energy_density(𝛙ₙ,j) for j ∈ 1:L], label = "|ψ$α⟩")
    end
    r = [parity(𝛙[α]) for α ∈ eachindex(𝛙)]
    cosk₀ = cosk(𝛙₀)

    𝛙ₛ = [(𝛙[α]+reflect(𝛙[α]))/(norm(𝛙[α]+reflect(𝛙[α]))) for α ∈ eachindex(𝛙)]
    rₛ = [parity(𝛙ₛ[α]) for α ∈ eachindex(𝛙)]
    println("Parity of the symmetrized states: ", rₛ)

    function trunccos(x)
        if x > 1
            return 1
        elseif x < -1
            return -1
        else
            return x
        end
    end

    k = [acos(trunccos(cosk(𝛙ₛ[α]))) for α ∈ eachindex(𝛙)]
    println("Momenta: ", k)

    # We select all the state with the parity which are the same of the first excited state
    # println("Number of states: ", length(𝛙))
    # 𝛙 = [𝛙[α] for α ∈ eachindex(𝛙) if parity(𝛙[α]) == parity(𝛙[1])]
    # println("Number of states with the same parity: ", length(𝛙))

    # We create a scatterplot of the energies vs momenta
    # scatter([k₀], [E₀], label = "E0")
    # scatter(k, E, label = "E(k)")
    # savefig("dispersion.png")

    scatter(k, E)
    savefig("dispersion.png")


    # Constructing the matrices
    H0 = MPO(opSumIsing(L, λ; n = 0, prefactor = 1/L), sites)
    H1 = MPO(opSumIsing(L, λ; n = 1, prefactor = 1/L), sites)
    H2 = MPO(opSumIsing(L, λ; n = 2, prefactor = 1/L), sites)
    A0 = inner(𝛙₀', H0, 𝛙₀)
    A1 = inner(𝛙₀', H1, 𝛙₀)
    A2 = inner(𝛙₀', H2, 𝛙₀)
    B0 = reduce(hcat, [[inner(𝛙[β]', H0, 𝛙[α]) for α ∈ eachindex(𝛙)] for β ∈ eachindex(𝛙)])
    B1 = reduce(hcat, [[inner(𝛙[β]', H1, 𝛙[α]) for α ∈ eachindex(𝛙)] for β ∈ eachindex(𝛙)])
    B2 = reduce(hcat, [[inner(𝛙[β]', H2, 𝛙[α]) for α ∈ eachindex(𝛙)] for β ∈ eachindex(𝛙)])

    # Assuming the matrix real, we find the eigenvalues and eigenvectors
    eig = eigen(B2)
    vals = eig.values
    vecs = eig.vectors
    zꜜ = vecs[:,1]
    println(vals[1])
    println("Coefficients: ", zꜜ)

    # Following steps if we cannot assume the matrix real
    #= function ℰ(n,θ)
        A = n==0 ? A0 : n==1 ? A1 : A2
        B = n==0 ? B0 : n==1 ? B1 : B2
        ℑ = eachindex(θ)
        return real(-L * A + sum(exp(im * (θ[α] - θ[β])) * B[α,β] for α ∈ ℑ, β ∈ ℑ))
    end
    σ²(θ) = ℰ(2,θ) / ℰ(0,θ) # - (ℰ(1,θ) / ℰ(0,θ))^2
    θ₀ = zeros(Float64, nₑ)
    result = optimize(σ², θ₀)
    θꜜ = result.minimizer
    println(result) =#

    # We sum all the 𝛙ₙ modulated by the z
    𝓌 = sum([zꜜ[α] * 𝛙[α] for α ∈ eachindex(𝛙)], cutoff = 10^(-15))

    # plot!([energy_density(𝓌,j) for j ∈ 1:L], label = "|W⟩")
    
    plot!([log(abs(energy_density(𝓌,j))) for j ∈ 1:L], label = "|W⟩")

    savefig("WannierLog.png")

    # We compute the outer product between the Wannier and the groundstate
    # A = outer(𝓌', 𝛙₀)
end