using ITensors

"""
Generates the Pauli matrices
"""
function pauli_matrix(symbol::String, 𝚒::Index, 𝚓::Index)
    pauli_matrix = ITensor(𝚒, 𝚓)

    if symbol == "x"
        pauli_matrix[𝚒=>1, 𝚓=>1] = 0 + 0im
        pauli_matrix[𝚒=>1, 𝚓=>2] = 1 + 0im
        pauli_matrix[𝚒=>2, 𝚓=>1] = 1 + 0im
        pauli_matrix[𝚒=>2, 𝚓=>2] = 0 + 0im
    elseif symbol == "y"
        pauli_matrix[𝚒=>1, 𝚓=>1] = 0 + 0im
        pauli_matrix[𝚒=>1, 𝚓=>2] = 0 - 1im
        pauli_matrix[𝚒=>2, 𝚓=>1] = 0 + 1im
        pauli_matrix[𝚒=>2, 𝚓=>2] = 0 + 0im
    elseif symbol == "z"
        pauli_matrix[𝚒=>1, 𝚓=>1] = 1 + 0im
        pauli_matrix[𝚒=>1, 𝚓=>2] = 0 + 0im
        pauli_matrix[𝚒=>2, 𝚓=>1] = 0 + 0im
        pauli_matrix[𝚒=>2, 𝚓=>2] = -1 + 0im
    elseif symbol == "+"
        pauli_matrix[𝚒=>1, 𝚓=>1] = 0 + 0im
        pauli_matrix[𝚒=>1, 𝚓=>2] = 1 + 0im
        pauli_matrix[𝚒=>2, 𝚓=>1] = 0 + 0im
        pauli_matrix[𝚒=>2, 𝚓=>2] = 0 + 0im
    elseif symbol == "-"
        pauli_matrix[𝚒=>1, 𝚓=>1] = 0 + 0im
        pauli_matrix[𝚒=>1, 𝚓=>2] = 0 + 0im
        pauli_matrix[𝚒=>2, 𝚓=>1] = 1 + 0im
        pauli_matrix[𝚒=>2, 𝚓=>2] = 0 + 0im
    else
        @warn "Symbol not valid for 𝛔; identity taken."
        pauli_matrix[𝚒=>1, 𝚓=>1] = 1 + 0im
        pauli_matrix[𝚒=>1, 𝚓=>2] = 0 + 0im
        pauli_matrix[𝚒=>2, 𝚓=>1] = 0 + 0im
        pauli_matrix[𝚒=>2, 𝚓=>2] = 1 + 0im
    end

    return pauli_matrix
end

function translation_operator(
    𝙸::Vector{Index{Int64}}, 
    𝙹::Vector{Index{Int64}}
    )::ITensor

    𝐓 = ITensor(𝙸, 𝙹)
    # Constructing the translation operator
    for ind_up in eachval(Tuple(x for x in 𝙸))
        arr = [x for x in Tuple(ind_up)]
        arr = circshift(arr, 1)
        tup = Tuple(arr)
        ind_down = CartesianIndex(tup)
        𝐓[CartesianIndex(ind_up, ind_down)] = 1
    end 
    return 𝐓
end

function ising_hamiltonian(
    𝙸::Vector{Index{Int64}}, 
    𝙹::Vector{Index{Int64}}, 
    J::Real, 
    h::Real
    )::ITensor

    𝐇 = ITensor(𝙸, 𝙹)
    𝛔ˣ(i) = pauli_matrix("x", 𝙸[i], 𝙹[i])
    𝛔ᶻ(i) = pauli_matrix("z", 𝙸[i], 𝙹[i])

    L = length(𝙸)
    for i in 1:(L-1)
        𝐇 += -h * 𝛔ˣ(i)
        𝐇 += -J/2 * 𝛔ᶻ(i) * 𝛔ᶻ(i+1)
    end
    𝐇 += -h * 𝛔ˣ(L)
    𝐇 += -J/2 * 𝛔ᶻ(L) * 𝛔ᶻ(1)
end

L = 3 # The number of sites
d = 2 # The local dimension

# We define the indices
𝙸 = siteinds(d, L)
𝙹 = siteinds(d, L)

Ci = combiner(𝙸, tags="ci")
Cj = combiner(𝙹, tags="cj")

𝐇 = ising_hamiltonian(𝙸, 𝙹, 1, 1)
𝐇′= Cj * 𝐇 * Ci

show(𝐇′)