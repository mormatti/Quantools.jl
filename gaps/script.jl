let
    d = 2
    L = 10
    X = [0 1; 1 0]
    Z = [1 0; 0 -1]
    I = [1 0; 0 1]
    M = multmatrix(d,L)

    function cutoff_entries(O,ϵ)
        O[abs.(O) .< ϵ] .= 0
        return O
    end

    function cutoff_operators(O, N)
        NOp = O
        # We set all the entries of O which corresponds to an entry of M which is > N to zero
        NOp[M .> N] .= 0
        return NOp
    end

    function gap_error(h; ϵ = 0.0, N = L)
        H = - Σⱼ(X,X,L) - h * Σⱼ(Z,L)
        F = eigen(H)

        gap_number = 2
        # n = λ > 0.5 ? gap_number+2 : gap_number+1
        n = gap_number + 1

        E0, E1 = F.values[1], F.values[n]
        ψ0, ψ1 = F.vectors[:,1], F.vectors[:,n]
        
        O = ψ1' ⊗ ψ0
        Δ = E1 - E0
        O = cutoff_entries(O,ϵ)
        O = cutoff_operators(O, N)
        Δ′ = compute_gap(H, O, ψ0, E0)
        return abs(Δ - Δ′) / abs(Δ)
    end

    λrange = 0:0.001:1.0
    logϵrange = -8:0.1:-1

    # We initialize a 2x2 matrix
    data = zeros(length(logϵrange), length(λrange))
    for i in eachindex(λrange)
        # We measure the time elalsed in the following loop
        δt = @elapsed for j in eachindex(logϵrange)
            data[j,i] = log(gap_error(λrange[i], ϵ = 10^logϵrange[j]))
        end
        println("Estimated remaining time: ", δt * (length(λrange) - i))
    end

    heatmap(λrange, logϵrange, data)
end