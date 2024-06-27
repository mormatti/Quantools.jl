"""
Exchanges the physical indices of a (finite) `MPS`.
This is possible only if the dimension of the physical indices is symmetric, i.e. 
if the dimension of a physical index at site `j` is the same as the dimension of the
physical index at site `N-j+1`, where `N` is the number of sites of the `MPS`.
"""
function reflect(ψ::MPS)
    N, Sd = length(ψ), siteinds(ψ)
    println(Sd)
    ϕ = MPS(N)
    for j in 1:N
        h = N - j + 1
        ϕ[h] = ψ[j] * delta(Sd[j], Sd[h]')
    end
    noprime!(siteinds, ϕ)
    return ϕ
end
export reflect