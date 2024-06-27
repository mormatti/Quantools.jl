mutable struct QuantumOperator
    parent_system   # QuantumSystem
    representation  # TensorialRepresentations
    time_instant    # Float64
    hermiticity     :: Bool
    unitarity       :: Bool
end

QuantumOperator(; quantum_system = missing, representation = missing, hermiticity = missing, unitarity = missing) = QuantumOperator(quantum_system, representation, hermiticity, unitarity)

parent_system(𝒪::QuantumOperator) = 𝒪.parent_system

isselfadjoint(𝒪::QuantumOperator) = isHermitian(𝒪)

isprojector(𝒪::QuantumOperator) = isIdempotent(𝒪) && selfAdjoint(𝒪)

mutable struct IdentityOperator
    parent_system   # QuantumSystem
    representation  # TensorialRepresentations
end

IdentityOperator(; quantum_system = missing, representation = missing) = IdentityOperator(quantum_system, representation)

isunitary(𝒪::IdentityOperator) = true

ishermitian(𝒪::IdentityOperator) = true

mutable struct HamiltonianOperator
    parent_system   # QuantumSystem
    representation  # TensorialRepresentations
    time_instant    # Float64
end

HamiltonianOperator(; quantum_system = missing, representation = missing) = HamiltonianOperator(quantum_system, representation)

isunitary(𝒪::HamiltonianOperator) = false

ishermitian(𝒪::HamiltonianOperator) = true

mutable struct EvolutionOperator
    parent_system   # QuantumSystem
    representation  # TensorialRepresentations
    time_instant     # Float64
end

isunitary(𝒪::EvolutionOperator) = true

ishermitian(𝒪::EvolutionOperator) = false

mutable struct MomentumOperator
    parent_system   # QuantumSystem
    representation  # TensorialRepresentations
end

MomentumOperator(; quantum_system = missing, representation = missing) = MomentumOperator(quantum_system, representation)

isunitary(𝒪::MomentumOperator) = false

ishermitian(𝒪::MomentumOperator) = true