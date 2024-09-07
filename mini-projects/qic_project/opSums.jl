# In this file we define the functions that return the Hamiltonian of the Hamiltonian
# of the ladder QED, written in the MPS formalism.

# we import the module ITensors
using ITensors

"This function returns an `OpSum` object which contains the Hamiltonian of the quantum 
transverse field Ising model. The constant J is the coupling constant between the spins, 
and the constant h represents the transverse magnetic field. N is the number of spins."
function opSumIsing(; N::Int64, J=1.0, h=0.0, λ=nothing)
    if λ !== nothing
        h = J * λ
    end
    os = OpSum()
    for j=1:N-1
        os += -J,"Sz",j,"Sz",j+1
    end
    for j=1:N
        os += -h,"Sx",j
    end
    return os
end

"""This function returns an `OpSum` object which contains the Hamiltonian of the QED on
a ladder. The constant g is the coupling constant between the photons and the spins, and
the constant L is the number of plaquettes of the ladder."""
function opSumLadder(; L, g=1.0)
    Nl = 3 * L + 1
    os = OpSum()

    # We add the electric field term to the Hamiltonian.
    for j=1:Nl
        os += g^2 / 2,"Sz",j,"Sz",j
    end
    # We add the magnetic field term to the Hamiltonian.
    for p=0:L-1
        j = 3 * p
        os += 1 / (2*g^2),"S-",j+1,"S-",j+2,"S+",j+3,"S+",j+4
        os += 1 / (2*g^2),"S+",j+1,"S+",j+2,"S-",j+3,"S-",j+4
    end
    return os
end

"""This function returns an `OpSum` object which contains the renormalized Hamiltonian 
of the QED on a ladder.
The constant λ is the analogue to the magnetic field in the Ising model, and
the constant L is the number of plaquettes of the ladder."""
function opSumLadderNormalized(; L, λ = 1.0)
    os = OpSum()
    # We add the electric field term to the Hamiltonian.
    for j=1:(3*L + 1)
        os += 2,"Sz",j,"Sz",j
    end
    # We add the magnetic field term to the Hamiltonian.
    for p=0:L-1
        j = 3 * p
        os += λ/2,"S-",j+1,"S-",j+2,"S+",j+3,"S+",j+4
        os += λ/2,"S+",j+1,"S+",j+2,"S-",j+3,"S-",j+4
    end
    return os
end

"""Same of the previous function, but with periodic boundary conditions."""
function opSumLadderNormalizedPBC(; L, λ = 1.0)
    os = OpSum()
    # We add the electric field term to the Hamiltonian.
    for j=1:(3*L)
        os += 2,"Sz",j,"Sz",j
    end
    # We add the magnetic field term to the Hamiltonian.
    for p=0:L-2
        j = 3 * p
        os += λ/2,"S-",j+1,"S-",j+2,"S+",j+3,"S+",j+4
        os += λ/2,"S+",j+1,"S+",j+2,"S-",j+3,"S-",j+4
    end
    j = 3 * (L-1)
    os += λ/2,"S-",j+1,"S-",j+2,"S+",j+3,"S+",1
    os += λ/2,"S+",j+1,"S+",j+2,"S-",j+3,"S-",1
    return os
end

function opSumMagneticTerm(; L)
    Nl = 3*L + 1
    os = OpSum()

    # We add the magnetic field term to the Hamiltonian.
    for p=0:L-1
        j = 3 * p
        os += "S-",j+1,"S-",j+2,"S+",j+3,"S+",j+4
        os += "S+",j+1,"S+",j+2,"S-",j+3,"S-",j+4
    end
    return os
end


"""This function returns an `OpSum` object which contains an operator which represents
the plaquette magnetic operator in the center of the ladder. The constant L is the
number of plaquettes of the ladder."""
function magneticOrderParameter(; L, λ=1.0)
    os = OpSum()
    P = Int(floor(L / 2))
    for p=0:P-1
        j = 3 * p
        os += λ/4,"S-",j+1,"S-",j+2,"S+",j+3,"S+",j+4
        os += λ/4,"S+",j+1,"S+",j+2,"S-",j+3,"S-",j+4
    end
    return os
end

"""This function returns an `OpSum` object which contains an operator which represents
the analogue of the magnetization in the Ising model. The constant L is the number of  
plaquettes of the ladder."""
function magnetization(; L)
    os = OpSum()
    P = Int(floor(L / 2))
    j = 3 * P
    k = 3 * (P + 1)
    os += "Sz",j+2,"Sz",k+2
    return os
end