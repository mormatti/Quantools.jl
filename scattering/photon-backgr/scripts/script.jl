let 
    # We define constants
    i1 = 1.0 + 0.0im
    𝟙 = [i1 0; 0 i1] # The identity 2x2 matrix
    𝛔ˣ = [0 i1; i1 0] # The first Pauli matrix
    𝛔ʸ = [0 -im; im 0] # The second Pauli matrix
    𝛔ᶻ = [i1 0; 0 -i1] # The third Pauli matrix
    𝛔⁺ = [0 i1; 0 0]
    𝛔⁻ = [0 0; i1 0]
    g⁴ = 1.0 # The fourth power of the coupling constant
    J = 1/2 # The coupling constant of the spin interaction
    hˣ = √8 / g⁴ # The transverse field
    hᶻ = 1 # The longitudinal field
    L = 9 # The number of sites of the chain
    N = 2^L # The number of states of the chain

    𝐇 = isingHamiltonian(J, hˣ, hᶻ, L)
    𝐓 = translationOperator(L, d)
    𝐡 = isingLocalHamiltonian(1, J, hˣ, hᶻ, L)
    (k, ℰ, 𝛙) = blochStates(𝐇, 𝐓)
    ℰ₀ = minimum(ℰ)
    𝛙₀ = 𝛙[argmin(ℰ)]
    k, ℰ, 𝛙 = firstBandStates(k, ℰ, 𝛙, L)
    𝓌 = find_wannier(𝛙, 𝐓, 𝐡, jᶜ, ℰ₀, L)
    plot(energyDensityArray(𝓌, 𝐡, 𝐓, L))
    savefig("data/energyWannier.png")
    𝐋 = [𝟙, 𝛔⁺]
    interpolate(𝛙₁, 𝛙₀, 𝐓, 3, 7, 𝐋, d, L)
end