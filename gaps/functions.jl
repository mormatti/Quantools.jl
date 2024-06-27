function compute_gap(H, O, ψ, E0)
    normalize(ψ)
    nϵ = O' * ψ
    oϵ = O * ψ
    nϵ2 = nϵ' * nϵ
    oϵ2 = oϵ' * oϵ
    Enϵ = nϵ' * H * nϵ
    Eoϵ = oϵ' * H * oϵ
    return (Enϵ + Eoϵ - E0 * (nϵ2 + oϵ2))/(nϵ2 - oϵ2)
end

function multmatrix(d,L)
    J::Matrix{Int64} = ones(d,d)
    K::Matrix{Int64} = J - Matrix(I, d, d)
    M::Matrix{Int64} = K
    for n in 2:L
        M = M ⊗ J + ones(d^(n-1),d^(n-1)) ⊗ K
    end
    return M
end