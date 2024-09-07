# This program computes the ground state energy per plaquette of the ladder
# with open boundary conditions, using the DMRG algorithm.

# We import the necessary packages
using ITensors
using LaTeXStrings

let 
    # We set the parameters of the simulation
    horizontalSpinRepresentation = "S=2"
    verticalSpinRepresentation = "S=4"
    L = 10 # The number of plaquettes of the ladder
    energyTolerance = 1E-5
    maxSweepsNumber = 1000000
    maxBondDimension = 300
    truncationErrorCutoff = 1E-6

    # We compute the number of links in the ladder
    Nl = 3 * L + 1 # The number of links in the ladder

    # We create the sites of the ladder
    sites = siteinds(horizontalSpinRepresentation, Nl)
    for i = 1:Nl
        if i % 3 == 1
            sites[i] = siteind(verticalSpinRepresentation)
        end
    end

    # We initialize the MPS as a random one
    ψin = randomMPS(sites, 150)

    # We create an array for the values of λ, from 0 to 3 with steps of 0.01
    λ_values = 0:(0.01/10):0.1

    # We create an array for the values of the ground state energy per plaquette
    energy_per_plaquette = zeros(length(λ_values))

    # We compute the ground state energy per plaquette for each value of λ
    for (i, λ) in enumerate(λ_values)
        os = opSumLadderNormalized(L = L, λ = λ)
        H = MPO(os, sites)
        observer = DMRGObserver(energy_tol = energyTolerance)

        # We compute the ground state energy per plaquette
        E0, ψ0 = dmrg(H, ψin;
            nsweeps = maxSweepsNumber, 
            maxdim = maxBondDimension, 
            cutoff = truncationErrorCutoff, 
            outputlevel = 1,
            observer = observer
            )
        E0N = E0 / L

        # We append to "data/DMRG/dsh1/L100/data.txt" the couple (λ, E0/L)
        exampleFileIOStream = open("data/DMRG/dsh4/L$L/bd300.txt","a")
        write(exampleFileIOStream, "$λ $E0N \n");
        close(exampleFileIOStream)
    end
end