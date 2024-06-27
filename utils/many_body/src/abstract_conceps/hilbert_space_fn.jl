# Constructors



# Getters

field(𝒮::HilbertSpace) = 𝒮.field

dimension(𝒮::HilbertSpace) = 𝒮.dimension

HilbertSpace(; field = missing, dimension = missing) = HilbertSpace(field, dimension)

name(𝒮::HilbertSpace) = "$(𝒮.field)^$(𝒮.dimension)"

ℂ(n) = HilbertSpace(field = "ℂ", dimension = n)

function ⊗(ℋ₁::HilbertSpace, ℋ₂::HilbertSpace)
    if ℋ₁.field != ℋ₂.field
        throw(ArgumentError("Impossible to perform the tensor product because the fields of the two Hilbert spaces are different."))
    end
    return HilbertSpace(field = ℋ₁.field, dimension = 𝒮₁.dimension * 𝒮₂.dimension)
end

qubit_hilbert_space = ℂ(2)

qutrit_hilbert_space = ℂ(2)

quditHilbertSpace(d) = ℂ(d)