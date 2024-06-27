module LinearAlgebraExtensions
    using LinearAlgebra

    include("diagonalization.jl")
    include("shortcuts.jl")
    include("matrix_construction.jl")
end

export LinearAlgebraExtensions