# A configuration with matter
mutable struct Config
    Lx :: Int
    Ly :: Int
    n  :: Int  # Zn symmetry
    sites  :: Matrix{Bool} # The matter sites, with index [x,y]
    xlinks :: Matrix{Int}  # The links -, with index [ix,y]
    ylinks :: Matrix{Int}  # The links |, with index [x,iy]
end
export Config

lengthx(C::Config) = C.Lx
export lengthx

lengthy(C::Config) = C.Ly
export lengthy

nparameter(C::Config)  = C.n
export nparameter

sites(C::Config)  = C.sites
export sites

xlinks(C::Config) = C.xlinks
export xlinks

ylinks(C::Config) = C.ylinks
export ylinks

n(C::Config) = C.n
export n

# We overload the == operator for the Config type
Base.:(==)(C1::Config, C2::Config) = (C1.Lx == C2.Lx && C1.Ly == C2.Ly && C1.n == C2.n && C1.sites == C2.sites && C1.xlinks == C2.xlinks && C1.ylinks == C2.ylinks)

function Config(Ly::Int, Lx::Int, n::Int) # m rows, n columns
    @assert Lx > 0 && Ly > 0
    sites  = falses(Lx, Ly)
    xlinks = zeros(Int8, Lx-1, Ly)
    ylinks = zeros(Int8, Lx, Ly-1)
    Config(Lx, Ly, n, sites, xlinks, ylinks)
end

function randomConfig(Ly::Int, Lx::Int, n::Int)
    @assert Lx > 0 && Ly > 0
    @assert n%2 == 1
    emax = (n-1)/2
    sites  = rand(0:1, Lx, Ly)
    xlinks = rand(-emax:emax, Lx-1, Ly)
    ylinks = rand(-emax:emax, Lx, Ly-1)
    Config(Lx, Ly, n, sites, xlinks, ylinks)
end
export randomConfig

function randomGaugeInvConfig(Ly::Int, Lx::Int, n::Int)
    @assert Lx > 0 && Ly > 0
    @assert n%2 == 1
    C = randomConfig(Lx, Ly, n)
    while !isGaugeInvariant(C)
        C = randomConfig(Lx, Ly, n)
    end
    return C
end
export randomGaugeInvConfig

function randomConfigWithoutMatter(Ly::Int, Lx::Int, n::Int)
    @assert Lx > 0 && Ly > 0
    @assert n%2 == 1
    emax = (n-1)/2
    sites = falses(Lx, Ly)
    for x in 1:Lx, y in 1:Ly
        if (x+y)%2 == 0
            sites[x, y] = true
        end
    end
    xlinks = rand(-emax:emax, Lx-1, Ly)
    ylinks = rand(-emax:emax, Lx, Ly-1)
    Config(Lx, Ly, n, sites, xlinks, ylinks)
end
export randomConfigWithoutMatter

function randomGaugeInvConfigWithoutMatter(Ly::Int, Lx::Int, n::Int)
    @assert Lx > 0 && Ly > 0
    @assert n%2 == 1
    C = randomConfigWithoutMatter(Lx, Ly, n)
    while !isGaugeInvariant(C)
        C = randomConfigWithoutMatter(Lx, Ly, n)
    end
    return C
end
export randomGaugeInvConfigWithoutMatter

function indomain(config::Config, x::Int, y::Int)
    Lx = config.Lx
    Ly = config.Ly
    return 1 ≤ y ≤ Ly && 1 ≤ x ≤ Lx
end

function indomain(config::Config, ix::Int, iy::Int, pos::String)
    # We control that the string is "L", "R", "T", "B"
    @assert pos == "L" || pos == "R" || pos == "T" || pos == "B"
    Lx = config.Lx
    Ly = config.Ly
    if pos == "L"
        return 1 ≤ ix ≤ Lx && 1 ≤ iy ≤ Ly-1
    elseif pos == "R"
        return 0 ≤ ix ≤ Lx-1 && 1 ≤ iy ≤ Ly-1
    elseif pos == "T"
        return 1 ≤ ix ≤ Lx-1 && 0 ≤ iy ≤ Ly-1
    elseif pos == "B"
        return 1 ≤ ix ≤ Lx-1 && 1 ≤ iy ≤ Ly
    end
end
export indomain

function Efield(config::Config, ix::Int, iy::Int, pos::String)
    if indomain(config, ix, iy, pos)
        if pos == "L"
            return config.ylinks[ix, iy]
        elseif pos == "R"
            return config.ylinks[ix+1, iy]
        elseif pos == "T"
            return config.xlinks[ix, iy+1]
        elseif pos == "B"
            return config.xlinks[ix, iy]
        end
    else
        return 0
    end
end
export Efield

function Eflux(config::Config, x::Int, y::Int) # Exiting electric flux
    @assert indomain(config, x, y)
    n = config.n
    E1 = Efield(config, x, y-1, "T")
    E2 = Efield(config, x, y, "L")
    E3 = Efield(config, x-1, y, "B")
    E4 = Efield(config, x-1, y-1, "R")
    return ZZ(n, E1 + E2 - E3 - E4)
end
export Eflux

function addEField(config::Config, ix::Int, iy::Int, pos::String, e::Int)
    n = config.n
    if indomain(config, ix, iy, pos)
        if pos == "L"
            config.ylinks[ix, iy] = ZZ(n, e + config.ylinks[ix, iy])
        elseif pos == "R"
            config.ylinks[ix+1, iy] = ZZ(n, e + config.ylinks[ix+1, iy])
        elseif pos == "T"
            config.xlinks[ix, iy+1] = ZZ(n, e + config.xlinks[ix, iy+1])
        elseif pos == "B"
            config.xlinks[ix, iy] = ZZ(n, e + config.xlinks[ix, iy])
        end
    end
end
export addEField

function isMatterSite(config::Config, x::Int, y::Int)::Bool
    @assert indomain(config, x, y)
    return (x+y)%2 == 0 ? true : false
end
export isMatterSite

function elecCharge(config::Config, x::Int, y::Int)
    @assert indomain(config, x, y)
    ms = isMatterSite(config, x, y)
    fd = config.sites[x, y]
    return ms ? (fd ? 0 : -1) : (fd ? +1 : 0)
end
export elecCharge

function isGaugeInvariant(config::Config)
    Lx = config.Lx
    Ly = config.Ly
    n = config.n
    for x in 1:Lx
        for y in 1:Ly
            if ZZ(n, elecCharge(config, x, y) - Eflux(config, x, y)) != 0
                return false
            end
        end
    end
    return true
end
export isGaugeInvariant

function withoutCharge(config::Config)
    Lx = config.Lx
    Ly = config.Ly
    for x in 1:Lx
        for y in 1:Ly
            if elecCharge(config, x, y) != 0
                return false
            end
        end
    end
    return true
end
export withoutCharge

function printZ3Config(config::Config)
    Lx = config.Lx
    Ly = config.Ly
    n = config.n
    @assert n == 3

    function symbx(e::Int)
        if e == 0
            return "-"
        elseif e == 1
            return "▶"
        elseif e == -1
            return "◀"
        end
    end

    function symby(e::Int)
        if e == 0
            return "|"
        elseif e == 1
            return "▲"
        elseif e == -1
            return "▼"
        end
    end

    function symbc(c::Int)
        if c == 0
            return "⋅"
        elseif c == 1
            return "⊕"
        elseif c == -1
            return "⊖"
        end
    end

    function printxline(y) # ○-○-○-○-○-○
        fsymbc(x,y) = symbc(elecCharge(config, x, y))
        fsymbx(x,y) = symbx(config.xlinks[x,y])
        for x in 1:Lx-1
            print(fsymbc(x,y), " ", fsymbx(x,y), " ")
        end
        println(fsymbc(Lx,y))
    end

    function printyline(y) # | | | | | |
        fsymby(x,y) = symby(config.ylinks[x,y])
        for x in 1:Lx-1
            print(fsymby(x,y), "   ")
        end
        println(fsymby(Lx,y))
    end

    for y in Ly:-1:1
        if y != Ly
            printyline(y)
        end
        printxline(y)
    end
end
export printZ3Config



# Operators

function opU(C::Config, ix::Int, iy::Int, m::Int) # U plaquette operator
    C′ = deepcopy(C)
    addEField(C′, ix, iy, "T", -m)
    addEField(C′, ix, iy, "L", -m)
    addEField(C′, ix, iy, "R", m)
    addEField(C′, ix, iy, "B", m)
    return 1.0, C′
end
export opU

function opψUψplusHc(C::Config, ix::Int, iy::Int, pos::String) # ψUψ† + H.c. operator
    C′ = deepcopy(C)
    if pos == "L"
        p1x, p1y, p2x, p2y = ix, iy, ix, iy+1
    elseif pos == "R"
        p1x, p1y, p2x, p2y = ix+1, iy, ix+1, iy+1
    elseif pos == "T"
        p1x, p1y, p2x, p2y = ix, iy+1, ix+1, iy+1
    elseif pos == "B"
        p1x, p1y, p2x, p2y = ix, iy, ix+1, iy
    end 
    c1 = elecCharge(C′, p1x, p1y)
    c2 = elecCharge(C′, p2x, p2y)
    C′.sites[p1x, p1y] = !C′.sites[p1x, p1y]
    C′.sites[p2x, p2y] = !C′.sites[p2x, p2y]
    c1′ = elecCharge(C′, p1x, p1y)
    c2′ = elecCharge(C′, p2x, p2y)
    Δ = (c2 - c2′) - (c1 - c1′)
    sign = Int(Δ / 2)
    addEField(C′, ix, iy, pos, Int(sign))
    coeff = abs(sign)
    return coeff, C′
end
export opψUψplusHc

function ψdagψ(C::Config, x::Int, y::Int) # ψ†ψ operator
    return abs(elecCharge(C, x, y)), C
end
export ψdagψ

function Esquared(C::Config, ix::Int, iy::Int, pos::String) # E² operator
    return Efield(C, ix, iy, pos)^2, C
end
export Esquared