# A configuration with matter
mutable struct ExtConfig
    physConfig :: Config
    vortexConfig :: Config
end
export ExtConfig

physConfig(C::ExtConfig) = C.physConfig
export physConfig

vortexConfig(C::ExtConfig) = C.vortexConfig
export vortexConfig

function lengthx(C::ExtConfig)
    Lx1 = lengthx(C.physConfig)
    Lx2 = lengthx(C.vortexConfig)
    @assert Lx1 == Lx2
    return Lx1
end
export lengthx

function lengthy(C::ExtConfig)
    Ly1 = lengthy(C.physConfig)
    Ly2 = lengthy(C.vortexConfig)
    @assert Ly1 == Ly2
    return Ly1
end
export lengthy

function n(C::ExtConfig)
    n1 = n(C.physConfig)
    n2 = n(C.vortexConfig)
    @assert n1 == n2
    return n1
end
export n

function Base.:(==)(C1::ExtConfig, C2::ExtConfig)
    C1.physConfig == C2.physConfig && C1.vortexConfig == C2.vortexConfig
end

function ExtConfig(Ly::Int, Lx::Int, n::Int) # m rows, n columns
    @assert Lx > 0 && Ly > 0
    ExtConfig(Config(Ly, Lx, n), Config(Ly, Lx, n))
end
export ExtConfig

function randomExtConfig(Ly::Int, Lx::Int, n::Int)
    @assert Lx > 0 && Ly > 0
    @assert n%2 == 1
    emax = (n-1)/2
    ExtConfig(randomConfig(Ly, Lx, n), randomConfigWithoutMatter(Ly, Lx, n))
end
export randomExtConfig

function randomGaugeInvExtConfig(Ly::Int, Lx::Int, n::Int)
    @assert Lx > 0 && Ly > 0
    @assert n%2 == 1
    return ExtConfig(randomGaugeInvConfig(Ly, Lx, n), randomGaugeInvConfigWithoutMatter(Ly, Lx, n))
end
export randomGaugeInvExtConfig

function isGaugeInvariant(C::ExtConfig)
    isGaugeInvariant(C.physConfig) && isGaugeInvariant(C.vortexConfig)
end
export isGaugeInvariant

function printZ3ExtConfig(C::ExtConfig)
    printZ3Config(C.physConfig)
    println("  ,")
    printZ3Config(C.vortexConfig)
end
export printZ3ExtConfig

function tensorProductBasis(Basis1, Basis2)
    Basis = []
    for C1 in Basis1, C2 in Basis2
        push!(Basis, ExtConfig(C1, C2))
    end
    return Basis
end
export tensorProductBasis

function mapping(C::ExtConfig)
    n = C.physConfig.n
    n′ = C.vortexConfig.n
    @assert n == n′

    C′ = deepcopy(C)

    Ax = C′.physConfig.xlinks
    Ay = C′.physConfig.ylinks
    Bx = C′.vortexConfig.xlinks
    By = C′.vortexConfig.ylinks

    for x in 1:size(Ax,1), y in 1:size(Ax,2)
        Ax[x, y] = ZZ(n, Ax[x, y] + Bx[x, y])
        Bx[x, y] = 0
    end
    for x in 1:size(Ay,1), y in 1:size(Ay,2)
        Ay[x, y] = ZZ(n, Ay[x, y] + By[x, y])
        By[x, y] = 0
    end

    return 1.0, C′
end
export mapping

function opU(C::ExtConfig, ix::Int, iy::Int, m::Int) # U plaquette operator
    C′ = deepcopy(C)
    r, C′.vortexConfig = opU(C.vortexConfig, ix, iy, m)
    return r, C′
end
export opU

function opψUψplusHc(C::ExtConfig, ix::Int, iy::Int, pos::String) # ψUψ† + H.c. operator
    C′ = deepcopy(C)
    r, C′.physConfig = opψUψplusHc(C.physConfig, ix, iy, pos)
    return r, C′
end
export opψUψplusHc

function ψdagψ(C::ExtConfig, x::Int, y::Int) # ψ†ψ operator
    return abs(elecCharge(C.physConfig, x, y)), C
end
export ψdagψ

function Esquared(C::ExtConfig, ix::Int, iy::Int, pos::String) # E² operator
    return (ZZ(n(C), Efield(C.physConfig, ix, iy, pos) + Efield(C.vortexConfig, ix, iy, pos)))^2, C
end
export Esquared

function HamiltonianMatrix(Basis, t::Real, m::Real, gE::Real, gB::Real)
    N = length(Basis)
    H = zeros(ComplexF64, N, N)
    Lx = lengthx(Basis[1])
    Ly = lengthy(Basis[1])
    for x ∈ 1:Lx, y ∈ 1:Ly # Mass term
        sign = (-1)^(x+y)
        H += sign * m * matrix(Basis, C -> ψdagψ(C, x, y))
    end
    for ix ∈ 1:Lx-1, iy ∈ 1:Ly-1 # Magnetic term
        H += -gB^2/2 * matrix(Basis, C -> opU(C, ix, iy, 1))
        H += -gB^2/2 * matrix(Basis, C -> opU(C, ix, iy, -1))
    end
    for ix ∈ 1:Lx-1, iy ∈ 1:Ly-1 # Electric term and hopping term
        H += -t * matrix(Basis, C -> opψUψplusHc(C, ix, iy, "L"))
        H += gE^2/2 * matrix(Basis, C -> Esquared(C, ix, iy, "L"))
        H += -t * matrix(Basis, C -> opψUψplusHc(C, ix, iy, "B"))
        H += gE^2/2 * matrix(Basis, C -> Esquared(C, ix, iy, "B"))
    end
    for ix ∈ 1:Lx-1
        iy = Ly-1
        H += -t * matrix(Basis, C -> opψUψplusHc(C, ix, iy, "T"))
        H += gE^2/2 * matrix(Basis, C -> Esquared(C, ix, iy, "T"))
    end
    for iy ∈ 1:Ly-1
        ix = Lx-1
        H += -t * matrix(Basis, C -> opψUψplusHc(C, ix, iy, "R"))
        H += gE^2/2 * matrix(Basis, C -> Esquared(C, ix, iy, "R"))
    end
    return H
end
export HamiltonianMatrix



