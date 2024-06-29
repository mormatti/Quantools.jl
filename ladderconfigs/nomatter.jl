using GraphRecipes, Plots
using LinearAlgebra

function printformathematica(matrix::Matrix)
    # Initialize an array to hold the string representation of each row
    rows = []
    
    # Iterate through each row of the matrix
    for row in eachrow(matrix)
        # Convert each element of the row to a string and join them with commas
        row_str = join(row, ", ")
        # Enclose the row string in curly braces and push to the rows array
        push!(rows, "{" * row_str * "}")
    end
    
    # Join all rows with commas and enclose the result in curly braces
    matrix_str = "{" * join(rows, ", ") * "}"
    
    # Print the final matrix string
    println(matrix_str)
end

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

function subtractZZ(n, a, b)
    return ZZ(n, a - b)
end

function subtractSpin(sx, sy, a, b)
    if a > sx || b > sx || a < -sx || b < -sx
        return nothing
    end
    c = a - b
    if c > sy || c < -sy
        return nothing
    end
    return c
end

"""Returns true if the plaquette B can follow the plaquette A"""
function follows(A, B, subtractfunction)
    XL = subtractfunction(B[1], A[1])
    XR = subtractfunction(A[2], B[2])
    if XL === nothing || XR === nothing
        return false
    end
    return XL == XR
end

function adMatrixZZ(N)
    max = (N-1)/2
    plaquettes = []
    for a in -max:1:max
        for b in -max:1:max
            plaquettes = push!(plaquettes, (a,b))
        end
    end
    M = zeros(Int, length(plaquettes), length(plaquettes))

    display(eachindex(plaquettes))

    for i ∈ eachindex(plaquettes)
        for j ∈ eachindex(plaquettes)
            subtractfunc(a,b) = subtractZZ(N, a, b)
            if follows(plaquettes[i], plaquettes[j], subtractfunc)
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
            plaquettes = push!(plaquettes, (a,b))
        end
    end

    display(plaquettes)

    M = zeros(Int, length(plaquettes), length(plaquettes))

    display(eachindex(plaquettes))

    for i ∈ eachindex(plaquettes)
        for j ∈ eachindex(plaquettes)
            subtractfunc(a,b) = subtractSpin(sx, sy, a, b)
            if follows(plaquettes[i], plaquettes[j], subtractfunc)
                M[i,j] = 1
            end
        end
    end

    return M
end




