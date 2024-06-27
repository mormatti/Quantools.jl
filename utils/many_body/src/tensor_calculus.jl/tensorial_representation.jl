using LinearAlgebra
using ITensors

mutable struct TensorialRepresentations
    vector         # Vector{ComplexF64}, StoredData
    matrix         # Matrix{ComplexF64}, StoredData
    tensor         # ITensor (or other), StoredData
    tensor_network # ITensor (or other), StoredData
end

TensorialRepresentations(; vector_form = missing, matrix_form = missing, tensor_form = missing, tensor_network_form = missing) = TensorialRepresentations(vector_form, matrix_form, tensor_form, tensor_network_form)

TensorialRepresentations(::Vector{ComplexF64}) = TensorialRepresentations(vector_form = vector_form, matrix_form = missing, tensor_form = missing, tensor_network_form = missing)

TensorialRepresentations(::Matrix{ComplexF64}) = TensorialRepresentations(vector_form = missing, matrix_form = matrix_form, tensor_form = missing, tensor_network_form = missing)

TensorialRepresentations(::ITensor) = TensorialRepresentations(vector_form = missing, matrix_form = missing, tensor_form = tensor_form, tensor_network_form = missing)

TensorialRepresentations(::MPS) = TensorialRepresentations(vector_form = missing, matrix_form = missing, tensor_form = missing, tensor_network_form = tensor_network_form)

function numberOfRepresentations(𝒯::TensorialRepresentations)::Int64
    n = 0
    if 𝒯.vector_form != missing
        n += 1
    end
    if 𝒯.matrix_form != missing
        n += 1
    end
    if 𝒯.tensor_form != missing
        n += 1
    end
    if 𝒯.tensor_network_form != missing
        n += 1
    end
    return n
end

isUnique(𝒯::TensorialRepresentations)::Bool = numberOfRepresentations(𝒯) == 1