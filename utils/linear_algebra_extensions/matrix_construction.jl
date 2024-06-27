# We create a function with the following inputs: an array A of type T; a function which 
# maps elements of type T in elements of type T preceded by a real number. The output is a
# matrix of complex numbers. 
function matrix(basis, func)
    N = length(basis)
    M = zeros(ComplexF64, N, N)
    for i in 1:N
        r, target = func(basis[i])
        if r != 0
            j = findfirst(item -> item == target, basis)
            M[j,i] = r
        end
    end
    return M
end

export matrix