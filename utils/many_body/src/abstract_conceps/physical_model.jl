mutable struct IsingModelParameters
    J   ::Float64
    hx  ::Float64
    hy  ::Float64
    hz  ::Float64
end

mutable struct PottsModelParameters
    J   ::Float64
    h   ::Float64
end

mutable struct HeisenbergModelParameters
    Jx   ::Float64
    Jy   ::Float64
    Jz   ::Float64
    hx   ::Float64
    hy   ::Float64
    hz   ::Float64
end

mutable struct HubbardModelParameters
    t   ::Float64
    U   ::Float64
    μ   ::Float64
end

abstract type IsingModel end

abstract type PottsModel end

abstract type HeisenbergModel end

SpinModel::Type = Union{IsingModel, PottsModel, HeisenbergModel}

abstract type BoseHubbardModel end

abstract type FermiHubbardModel end

HubbardModel::Type = Union{BoseHubbardModel, FermiHubbardModel}

DiscreteModel::Type = Union{IsingModel, PottsModel, HeisenbergModel, HubbardModel}

abstract type ContinuousModel end

PhysicalModel::Type = Union{DiscreteModel, ContinuousModel}