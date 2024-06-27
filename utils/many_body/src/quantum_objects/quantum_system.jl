mutable struct ElementarySystem
    parent_system    # QuantumSystem
    hilbert_space    # HilbertSpace
    properties       # Vector{Properties} or Any
end

space_dimension(𝒮::pointSystem) = 0

mutable struct CompositeSystem
    parent_system               # QuantumSystem
    hilbert_space               # HilbertSpace
    subsystems                  # Vector{QuantumSystems}
    properties                  # Vector{Properties} or Any

    real_space                  # RealSpace
    physical_model              # PhysicalModel
end

ManyBodySystem::Type = CompositeSystem

ManyBodySystem(; name = missing, physical_model = missing, parameters = missing, hilbert_space = missing, physical_space = missing, subsystems = missing, boundary_conditions = missing, dynamical_properties = missing, translational_properties = missing) = ManyBodySystem(name, physical_model, parameters, hilbert_space, physical_space, subsystems, boundary_conditions, dynamical_properties, translational_properties)

space_dimension(𝒮::ManyBodySystem) = spaceDimension(𝒮.real_space)

mutable struct SpacetimeProperties
    parent_system   # QuantumSystem
    spacetime       # Spacetime
    time_operator   # QuantumOperator
    time_eigenbasis # Array{QuantumStates}
end

mutable struct DynamicalProperties
    parent_system       # QuantumSystem
    hamiltonian         # QuantumOperator
    boundary_conditions # BoundaryConditions
    evolution_operator  # QuantumOperator
    energy_eigenbasis   # Array{QuantumStates}
    groundstate_basis   # QuantumState
end
   
DynamicalProperties(; quantum_system = missing, hamiltonian = missing, evolution_operator = missing, energy_eigenbasis = missing, groundstate_basis = missing) = DynamicalProperties(quantum_system, hamiltonian, evolution_operator, energy_eigenbasis, groundstate_basis)

mutable struct TranslationalProperties
    parent_system           # QuantumSystem
    reciprocal_space        # ReciprocalSpace
    momentum_operator       # QuantumOperator
    translation_operator    # QuantumOperator
    momentum_eigenbasis     # Array{QuantumStates}
    period                  # Int64
    direction               # Int64
end

TranslationalProperties(; quantum_system = missing, momentum_operator = missing, translation_operator = missing, momentum_eigenbasis = missing, period = missing, direction = missing) = TranslationalProperties(quantum_system, momentum_operator, translation_operator, momentum_eigenbasis, period, direction)

QuantumSystem = Union{pointSystem, ManyBodySystem}