abstract type DirichletBoundaryConditions end

DBC::Type = DirichletBoundaryConditions

abstract type NeumannBoundaryConditions end

NBC::Type = NeumannBoundaryConditions

OpenBoundaryConditions = Union{DBC, NBC}

OBC::Type = OpenBoundaryConditions

abstract type PeriodicBoundaryConditions end

PBC::Type = PeriodicBoundaryConditions

BoundaryConditions::Type = Union{OBC, PBC}

