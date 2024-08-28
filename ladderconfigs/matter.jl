using GraphRecipes, Plots
using LinearAlgebra

"""A shortcut binary notation for the periodic modulus. Takes two integers."""
function ↻(n, m)
    n > 0 ? (n-1)%m + 1 : m + n%m
end

"""A shortcut binary notation for the periodic modulus centered in zero. Takes two integers."""
function ZZ(n, p)
    s = (n+1)/2
    p = p + s
    return (p > 0 ? (p-1)%n + 1 : n + p%n) - s
end

function is_same_kind(a, b)
    # Check if both are integers
    if isInteger(a) && isInteger(b)
        return true
    end

    # Check if both are semi-integers
    if (a % 1 == 0.5) && (b % 1 == 0.5)
        return true
    end

    # If neither condition is met, they are not of the same kind
    return false
end

# Helper function to check if a number is an integer
function isInteger(x)
    return x % 1 == 0
end

function spinSelect(s, p)
    if p > s || p < -s || !is_same_kind(s, p)
        return nothing
    end
    return p
end

"""Returns true if the plaquette B can follow the plaquette A"""
function follows(A, B, selection)
    X = selection(-A[1] + B[2] - A[3])
    if X === nothing
        return false
    else
        charge = selection(-A[2] + B[1] + X)
        if charge ∈ [-1,0]
            return true
        else
            return false
        end
    end
end

function adMatrixZZ(N)
    max = (N-1)/2
    plaquettes = []
    for a in -max:1:max
        for b in -max:1:max
            for c in [0,1]
                plaquettes = push!(plaquettes, (a,b,c))
            end
        end
    end
    M = zeros(Int, length(plaquettes), length(plaquettes))

    display(eachindex(plaquettes))

    for i ∈ eachindex(plaquettes)
        for j ∈ eachindex(plaquettes)
            selection(a) = ZZ(N, a)
            if follows(plaquettes[i], plaquettes[j], selection)
                M[i,j] = 1
            end
        end
    end

    return M
end

function adMatrixSpin(sx, sy)
    plaquettes = []
    for a in -sx:1:sx
        for b in -sx:1:sx
            for c in [0,1]
                plaquettes = push!(plaquettes, (a,b,c))
            end
        end
    end

    display(plaquettes)

    M = zeros(Int, length(plaquettes), length(plaquettes))

    display(eachindex(plaquettes))

    for i ∈ eachindex(plaquettes)
        for j ∈ eachindex(plaquettes)
            selection(a) = spinSelect(sy, a)
            if follows(plaquettes[i], plaquettes[j], selection)
                M[i,j] = 1
            end
        end
    end

    return M
end




