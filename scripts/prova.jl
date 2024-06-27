using ITensors

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

let 
    L = 21
    λ = 0.5
    prefactor = 1.0
    H1 = opSumIsing(L, λ)
    H2 = opSumIsingPBC(L, λ)
    sites = siteinds("S=1/2", L)
    H1 = MPO(H1, sites)
    H2 = MPO(H2, sites)

    ψ = randomMPS(sites, 5)
    println(linkdims(ψ))
    println(linkdims(translate(translate(ψ))))
end