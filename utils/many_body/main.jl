function IdentityOperator(𝒮::QuantumSystem)
    D = 𝒮.hilbert_space.dimension
    𝐈 = Matrix{ComplexF64}(I, D, D)
    return IdentityOperator(quantum_system = 𝒮, representation = 𝐈)
end

random_state(𝒮::QuantumSystem) = QuantumState(quantum_system = 𝒮, representation = random_state_representation(𝒮), missing)