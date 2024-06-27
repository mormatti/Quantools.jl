module QuantumProjects
    # Utils
    include("utils/main.jl")

    # Projects
    include("gaps/main.jl")
    include("magnetic/main.jl")
    include("scattering/main.jl")
end

export QuantumProjects