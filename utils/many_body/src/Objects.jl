using LinearAlgebra
using ITensors

mutable struct HilbertSpace
    dimension::Int
end

mutable struct QuantumSystem
    # Properties
    hilbert_space::HilbertSpace
    hamiltonian_operator::Matrix{Complex}
end

mutable struct ManyBodySystem
    # Inherit from
    quantum_system::QuantumSystem
    # Properties
    number_of_local_systems::Int
    physical_indices::Vector{Index}
end

mutable struct IsingModel
    # Inherit from
    many_body_system::ManyBodySystem
    # Properties
    coupling_constant::Real
    longitudinal_field::Real
    transverse_field::Real
end