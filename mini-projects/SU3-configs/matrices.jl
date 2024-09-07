using LinearAlgebra

module QCD2

    using LinearAlgebra

    function printForMathematica(M)
        print("{")
        for i in 1:size(M, 1)
            print("{")
            for j in 1:size(M, 2)
                print(M[i, j])
                if j < size(M, 2)
                    print(", ")
                end
            end
            print("}")
            if i < size(M, 1)
                print(", ")
            end
        end
        print("}")
    end

    diagMpure = [0, 2, 2, 1, 2, 1, 1, 2, 1,
            1, 1, 2, 1, 1, 2, 2, 1, 1,
            2, 0, 2, 2, 0, 2, 2, 2, 2,
            0, 2, 2, 2, 2, 0, 2, 2, 1,
            2, 1, 1, 2, 2, 1, 1, 2, 1,
            1, 1, 2, 1, 1, 2, 0, 2, 2]
    
    Mpure = diagm(diagMpure)

    diagMups = [0, 0, 0, 1, 1, 1, 0, 0, 0,
            2, 2, 2, 1, 1, 1, 1, 0, 0, 
            0, 3, 3, 3, 2, 2, 2, 2, 2,
            1, 1, 1, 1, 1, 0, 0, 0, 3,
            3, 3, 2, 2, 2, 2, 1, 1, 1, 
            3, 3, 3, 2, 2, 2, 3, 3, 3]
    
    Mups = diagm(diagMups)

    function generateAudagger()
        Aud = zeros(Float64, 54, 54)
        Aud[4,1] = +√3
        Aud[5,2] = -√2
        Aud[6,3] = +1
        Aud[10,4] = +2
        Aud[11,5] = -√2
        Aud[12,6] = +√2
        Aud[13,7] = +√2
        Aud[14,8] = -1
        Aud[15,9] = +1
        Aud[16,9] = -√2
        Aud[20,10] = +√3
        Aud[21,11] = +1
        Aud[22,12] = +√2
        Aud[23,13] = +√2
        Aud[24,14] = +√(2/3)
        Aud[25,14] = (-√2)*√(2/3)
        Aud[26,15] = (√2)*√(2/3)
        Aud[27,15] = (√2)*√(2/3)
        Aud[27,16] = -√3
    end

    Aud = generateAudagger()

    # Matrices for Z3 link symmetry

    η = exp(2*im*π/3)

    diagWL = [1, η', η, η', η, 1, η', η, 1,
            η, 1, η', η, 1, η', η', η, 1,
            η', 1, η', η, 1, η', η', η, η,
            1, η', η', η, η, 1, η', η, η',
            η, 1, η', η, η, 1, η', η, 1,
            η, 1, η', η, 1, η', 1, η', η]

    WL = diagm(diagWL)

    diagWR = [1, η, η', 1, η, η', 1, η, η',
            1, η, η', 1, η, η', η', 1, η,
            η', 1, η, η', 1, η, η, η', η',
            1, η, η, η', η', 1, η, η', 1,
            η, η', 1, η, η, η', 1, η, η',
            1, η, η', 1, η, η', 1, η, η']

    WR = diagm(diagWR)

    function generateAdjecencyMatrix()
        N = length(WL)
        global M = zeros(Integer, N, N)
        for i in eachindex(WL), j in eachindex(WR)
            if WL[i] * WR[j] ≈ 1
                M[i, j] = 1
            end
        end
        return M
    end

    Adjac = generateAdjecencyMatrix()
end

let 
    QCD2.printForMathematica(QCD2.Mpure)
end