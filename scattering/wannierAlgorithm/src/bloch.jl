


"""
Finds the bloch states of the first band.

Inputs:
- `k` is the `Vector` of momenta of the eigenstates;
- `ℰ` is the list of energies of the eigenstates;
- `𝛙` is the `Matrix` of eigenvectors of the eigenstates.
- `L` is the number of sites of the chain.

Assumptions:
- The groundstate must be non-degenerate.
"""
function first_band_states(
    k::Vector{Float64},
    ℰ::Vector{Float64},
    𝛙::Vector{Vector{ComplexF64}},
    L::Integer
    )::Tuple{Vector{Float64}, Vector{Float64}, Vector{Vector{ComplexF64}}}

    @debug "Computing the first band states..."

    # Constructing the vector of tuples (k, ℰ, 𝛙)
    𝓣 = [(k[i], ℰ[i], 𝛙[i]) for i in eachindex(k)]

    # Selecting the groundstate and deleting it
    j = argmin([t[2] for t in 𝓣]) # the groundstate index for 𝓣
    deleteat!(𝓣, j) # Deleting the element of states at that index

    # Selecting the first band states
    # The list of the tuple states of the first band
    𝓣′ = []
    while length(𝓣) > 0
        # The index and momentum of the state with the minimum energy
        j′ = argmin([t[2] for t in 𝓣])
        k′ = 𝓣[j′][1]
        # Filling 𝓣′ with the state with the minimum energy
        push!(𝓣′, 𝓣[j′])
        # deleting all the elements with same (similar) momentum
        𝓣 = [t for t in 𝓣 if abs(k′ - t[1]) > π/L]
        # we repeat the process until 𝓣 is empty
    end

    # Transforming the a vector of tuples 𝓣′ into a tuple of vectors 𝓑′
    # the vector of momenta of the first band
    k = [t[1] for t in 𝓣′]
    # the vector of energies of the first band
    ℰ = [t[2] for t in 𝓣′]
    # the matrix of eigenvectors of the first band
    𝛙 = [t[3] for t in 𝓣′]

    # Returning the tuple of bloch states of the first band
    return k, ℰ, 𝛙
end
export firstBandStates

""" Plots the energy density profile of the state `𝛙`.

Inputs:
- `𝛙` is a generic qualtum state;
- `𝐡₁` is the local Hamiltonian of the first site;
- `𝐓` is the translation operator of the system;
- `L` is the number of sites of the chain.

Optional inputs:
- `ℰ₀` is the groundstate energy of the system (default 0.0).

Outputs:
- The `Vector{Real}` of the energy density profile of the state `𝛙`.
"""
function energy_density_array(
    𝛙::Vector{ComplexF64}, 
    𝐡₁::Matrix{ComplexF64}, 
    𝐓::Matrix{ComplexF64}, 
    L::Integer;
    ℰ₀::Float64 = 0.0,
    )::Vector{Float64}

    @debug "Computing the energy density profile..."

    # We initialize the array
    arr::Vector{Float64} = [real(𝛙' * 𝐡₁ * 𝛙 - ℰ₀/L)]

    # The loop to compute all the others energies
    for _ in 2:L
        𝛙 = 𝐓' * 𝛙
        push!(arr, real(𝛙' * 𝐡₁ * 𝛙 - ℰ₀/L))
    end

    return arr
end
export energy_density_array