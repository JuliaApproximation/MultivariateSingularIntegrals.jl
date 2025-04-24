using MultivariateSingularIntegrals, QuadGK, ClassicalOrthogonalPolynomials, Test


Z = (2, 2im, 2+2im, 2-2im, -2 - 2im, -1.001-1.5im, -0.999-1.5im, -0.5-2im, -1-1.5im, -2+2im, -2, 0.1+0.2im, 0.1+im, 0.1-im, 1+0.1im, -1+0.1im, 1+im, -1+im, -1-im, 1-im)
n = 6
for z in Z
    L = logkernelsquare(z,n)
    for j = 0:n-1, k=0:n-j-1
        q = quadgk(s -> quadgk(t -> iszero(s+im*t) ? 0.0+0.0im : log(z-(s+im*t)) * legendrep(k,s) * legendrep(j,t), -1, 1, atol=1E-12)[1], -1, 1, atol=1E-12)[1]
        @test q ≈ L[k+1,j+1] atol=1E-10
    end
end

