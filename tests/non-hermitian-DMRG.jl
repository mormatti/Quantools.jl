using ITensors

let 
    sites = siteinds("S=1/2", 5)
    ψ₀ = randomMPS(sites, 3)

    H = randomMPO(sites)

    E₀, ψ₀ = dmrg(H, randomMPS(sites, 3); ishermitian=false, nsweeps=10)
    E₁, ψ₁ = dmrg(H, [ψ₀], randomMPS(sites, 3); ishermitian=false, nsweeps=10)

    println("E₀ = ", E₀)
    println("E₁ = ", E₁)
    println("|E₀| = ", abs(E₀))
    println("|E₁| = ", abs(E₁))
end