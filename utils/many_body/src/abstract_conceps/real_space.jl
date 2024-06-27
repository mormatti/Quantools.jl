mutable struct SetOfPoints
    points
end

mutable struct OpenLinearChain
    separation      ::Float64
    number_of_sites ::Int64
end

cardinality(𝒞::OpenLinearChain) = 𝒞.number_of_sites

mutable struct ClosedLinearChain
    separation      ::Float64
    number_of_sites ::Int64
end

cardinality(𝒞::ClosedLinearChain) = 𝒞.number_of_sites

mutable struct InfiniteLinearChain
    separation      ::Float64
end

cardinality(::InfiniteLinearChain) = Inf

LinearChain = Union{OpenLinearChain, ClosedLinearChain}

abstract type SquareLattice end

abstract type HexagonalLattice end

abstract type RectangularLattice end

abstract type CubicLattice end

isLinearChain(::LinearChain) = true

Lattice1D = Union{LinearLattice}

dimension(::Lattice1D) = 1

Lattice2D = Union{HexagonalLattice, RectangularLattice, SquareLattice}

dimension(::Lattice2D) = 2

Lattice3D = Union{CubicLattice}

dimension(::Lattice3D) = 3

Lattice = Union{Lattice1D, Lattice2D, Lattice3D}

DiscreteSpace = Union{SetOfPoints, Lattice}

#---------------------------------------------------------------

mutable struct InfiniteStraightLine end

mutable struct Circle
    length ::Float64
end

mutable struct InfinitePlane end

mutable struct FlatTorus
    length_x ::Float64
    length_y ::Float64
end

ContinuumSpace = Union{SetOfPoints, Lattice}

#---------------------------------------------------------------

RealSpace = {DiscreteSpace, ContinuumSpace}