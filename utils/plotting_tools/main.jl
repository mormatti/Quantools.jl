module Plotting
    using Plots

    """
    Plots a dispersion relation-like scatter-plot.

    Inputs:
    - `k` is the `Vector` of momenta of the eigenstates;
    - `ℰ` is the list of energies of the eigenstates;
    - `filename` is the name of the file (with extension) where the plot will be saved;
    - `aspectRatio` is the ratio between the width and the height of the plot (default 1:1).
    """
    function dispersion_relation(
        k::Vector{Float64},
        ℰ::Vector{Float64},
        filename::String = "plot.png", 
        aspectRatio::Float64 = 1.0)

        @debug "Plotting the dispersion relation..."

        # Plotting the dispersion relation
        plot(k, ℰ, seriestype=:scatter, markersize=4, legend=false, xlabel="k", ylabel="ℰ")
        plot!(size=(170,170 * aspectRatio), dpi=300) # Setting the size of the plot
        savefig(filename) # Saving the plot in a file
    end
    export dispersion_relation

    """
    Plots a matrix `M` as a heatmap.

    Inputs:
    - `M` is the matrix to be plotted;
    - `filename` is the name of the file (with extension) where the plot will be saved;
    - `dpi` is the resolution of the plot (default 300).
    """
    function matrix(M::Matrix; filepath::String = "matrixPlot.png", dpi::Integer = 300)
        M = abs.(M) # we convert the matrix elements to their modulus
        M = M[end:-1:1,:] # We reverse all the y of M 
        plot(heatmap(M, aspect_ratio=:equal), dpi = dpi)
        savefig(filepath)
    end
    export matrix

end

export Plotting