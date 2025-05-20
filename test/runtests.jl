using MultivariateSingularIntegrals, QuadGK, ClassicalOrthogonalPolynomials, LinearAlgebra, DoubleFloats, Test
import MultivariateSingularIntegrals: M1

@testset "L and M" begin
    @test M1(-2im) ≈ M1(-2im+eps()) ≈ M1(-2im-eps())
end

Z̃ = (2.0, 2.0im, 2.0+2im, 2.0-2im, -2.0 - 2im)
Z = (2, 2im, 2+2im, 2-2im, -2 - 2im, -1.001-1.5im, -0.999-1.5im, -0.5-2im, -1-1.5im, -2+2im, -2, 0.1+0.2im, 0.1+im, 0.1-im, 1+0.1im, -1+0.1im, 1+im, -1+im, -1-im, 1-im)
@testset "stieltjes" begin
    @testset "accurate for low order" begin
        n = 6
        for z in Z̃
            L = stieltjessquare(z,n)
            for j = 0:n-1, k=0:n-j-1
                q = quadgk(s -> quadgk(t -> iszero(s+im*t) ? 0.0+0.0im : inv(z-(s+im*t)) * legendrep(k,s) * legendrep(j,t), -1, 1, atol=1E-12)[1], -1, 1, atol=1E-12)[1]
                @test q ≈ L[k+1,j+1] atol=1E-10
            end
        end
    end

    @testset "BigFloat" begin
        setprecision(2000) do
            for z in big.(Z̃)
                P = Legendre{BigFloat}()
                n = 100
                M = Diagonal((P'P)[1:n,1:n])
                @test P[big(0.), 1:n]' *  inv(M)  * stieltjessquare(z, n) * inv(M) * P[big(0.), 1:n] ≈ inv(z) atol=1E-40
            end
        end
    end
end

@testset "logkernel" begin
    @testset "accurate for low order" begin
        n = 6
        for z in Z
            L = logkernelsquare(z,n)
            for j = 0:n-1, k=0:n-j-1
                q = quadgk(s -> quadgk(t -> iszero(s+im*t) ? 0.0+0.0im : log(z-(s+im*t)) * legendrep(k,s) * legendrep(j,t), -1, 1, atol=1E-12)[1], -1, 1, atol=1E-12)[1]
                @test q ≈ L[k+1,j+1] atol=1E-10
            end
        end
    end

    @testset "BigFloat" begin
        setprecision(2000) do
            for z in big.(Z̃)
                P = Legendre{BigFloat}()
                n = 100
                M = Diagonal((P'P)[1:n,1:n])
                @test P[big(0.), 1:n]' *  inv(M)  * logkernelsquare(z, n) * inv(M) * P[big(0.), 1:n] ≈ log(z) atol=1E-40
            end

            @test logkernelsquare(Double64(0) + 0im, 10) ≈ logkernelsquare(big(0.0) + 0im, 10) atol=1E-30
        end
    end
end

# n = 4000
# z = 2.0
# @profview logkernelsquare(z, n)

# x = range(-1,1,n); @time log.(abs.(x .+ im .* x'));