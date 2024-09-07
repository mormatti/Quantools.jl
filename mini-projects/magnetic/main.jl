module Magnetic
    # using
    using LinearAlgebra
    include("../utils/main.jl")
    using .Utils
    using .BasicTools
    using .LinearAlgebraExtensions

    # src
    include("src/Configs.jl")
    include("src/ExtConfigs.jl")

    # scripts
    # include("scripts/script1.jl")
    # include("scripts/script2.jl")
end

export Magnetic