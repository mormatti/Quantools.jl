let
    T = translationOperator(5)
    H1 = isingLocalHamiltonian(1, 1, 1, 1, 5)
    H2 = isingLocalHamiltonian(2, 1, 1, 1, 5)
    println(norm(T' * H1 * T - H2))
end