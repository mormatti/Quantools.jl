mutable struct QuantumState
    quantum_system   # QuantumSystem
    representation   # TensorialRepresentations
    time_instant     # Float64
    symmetries       # Array{Symmetry}
end

QuantumState(; quantum_system = missing, representation = missing, symmetries = missing) = QuantumState(quantum_system, representation, time_instant, symmetries)