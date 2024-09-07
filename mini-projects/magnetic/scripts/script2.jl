let 
    println("_____________")

    C = randomConfig(2, 2, 3)
    found = false
    while !found
        C = randomConfig(2, 2, 3)
        if isGaugeInvariant(C)
            found = true
        end
    end

    printZ3Config(C)
    println()

    C′ = opψUψplusHc(C, 1, 1, "B")[2]
    printZ3Config(C′)
    println()
end

let
    St = [true, false]
    Lk = [-1,0,1]
    Basis = []

    # We do the same loop of before but considering a 3x2 lattice
    for s1 in St, s2 in St, s3 in St, s4 in St, l1 in Lk, l2 in Lk, l3 in Lk, l4 in Lk
        C = Config(2, 2, 3)
        C.sites[1,1] = s1
        C.sites[1,2] = s2
        C.sites[2,1] = s3
        C.sites[2,2] = s4
        C.xlinks[1,1] = l1
        C.xlinks[1,2] = l2
        C.ylinks[1,1] = l3
        C.ylinks[2,1] = l4
        if isGaugeInvariant(C)
            push!(Basis, C)
        end
    end

    println("Number of configurations: ", length(Lst))
    N = length(Lst)

    # # We print the first 10 gauge invariant configurations
    # for i in 1:min(100, length(Lst))
    #     printZ3Config(Lst[i])
    #     println()
    # end

    H = zeros(ComplexF64, N, N)

    matrixPlot(Matrix(Lst, C -> Esquared(C, 1, 1, "R")))
end