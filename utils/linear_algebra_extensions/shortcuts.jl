# Kronecker product
function ⊗(a,b)
    return kron(a,b)
end

export ⊗

function 𝕀(d,n)
    return Matrix(I, d^n, d^n)
end

export 𝕀

function Σⱼ(aⱼ,L)
    d = size(aⱼ,1)
    𝕀d = 𝕀(d,1)
    Pₙ = aⱼ
    for n in 2:L
        Pₙ = Pₙ ⊗ 𝕀d + 𝕀(d,n-1) ⊗ aⱼ
    end
    return Pₙ
end

function Σⱼ(aⱼ,bⱼ₊₁,L)
    d = size(aⱼ,1)
    𝕀d = 𝕀(d,1)
    Pₙ = aⱼ ⊗ bⱼ₊₁
    for n in 3:L
        Pₙ = Pₙ ⊗ 𝕀d + 𝕀(d,n-2) ⊗ (aⱼ ⊗ bⱼ₊₁)
    end
    return Pₙ
end

export Σⱼ