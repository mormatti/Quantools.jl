let
    for g4M in [0.01, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
        using Plots
        plot(phases, moduli, seriestype=:scatter, markersize=4, legend=false, xlabel="Momentum", ylabel="Energy")
        plot!(size=(170,1000), dpi=300)
        fontsize = 12
        plot!(xticksfontsize=fontsize, yticksfontsize=fontsize, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontsize=fontsize)
        title!("g^4 = $g4M")
        
        # We save the plot in a file
        savefig("fileg4=$g4M.png")
    end
        
        
    # We concatenate all the pictures "fileg4=$g4M.png" in a single picture
    using Images
    using FileIO
    using ImageMagick
        
    # We create an array of images
    images = []
    for g4M in [0.01, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
        push!(images, load("fileg4=$g4M.png"))
    end
        
    # We concatenate the images
    image = hcat(images...)
        
    # We save the image
    save("file.png", image)
end