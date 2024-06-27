let # Beginning of script
     using ArnoldiMethod, LinearAlgebra, SparseArrays
     A = spdiagm(
          -3 => fill(-5.0, 20000),
         -1 => fill(1.0, 6000),
          0 => fill(-2.0, 700), 
          1 => fill(1.0, 6000),
          3 => fill(-5.0, 20000)
     )
     display(A)
     decomp, history = partialschur(A, nev=100, tol=1e-10, which=SR())
     display(history)
     λs, X = partialeigen(decomp)
     display(λs)
# end of script
end