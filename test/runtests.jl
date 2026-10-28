using MultivariateSingularIntegrals, QuadGK, ClassicalOrthogonalPolynomials, LinearAlgebra, DoubleFloats, Test
import MultivariateSingularIntegrals: M1

@testset "L and M" begin
    @test M1(-2im) ≈ M1(-2im+eps()) ≈ M1(-2im-eps())
end

Z̃ = (2.0, 2.0im, 2.0+2im, 2.0-2im, -2.0 - 2im)
Z = (2, 2im, 2+2im, 2-2im, -2 - 2im, -1.001-1.5im, -0.999-1.5im, -0.5-2im, -1-1.5im, -2+2im, -2, 0.1+0.2im, 0.1+im, 0.1-im, 1+0.1im, -1+0.1im, 1+im, -1+im, -1-im, 1-im, 1-0.3im, -1-0.3im, 1.0-0.0im, -1.0-0.0im)
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

@testset "logkernel square" begin
    @testset "accurate for low order" begin
        n = 6
        for z in Z
            L = logkernelrectangle(z,n)
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
                @test P[big(0.), 1:n]' *  inv(M)  * logkernelrectangle(z, n) * inv(M) * P[big(0.), 1:n] ≈ log(z) atol=1E-40
            end

            @test logkernelrectangle(Double64(0) + 0im, 10) ≈ logkernelrectangle(big(0.0) + 0im, 10) atol=1E-30
        end
    end
end

@testset "logkernelrectangle" begin
    @testset "accurate for low order" begin
        n = 6
        for γ in (1.3, 0.6)
            for z in (2γ, 2im, 2γ+2im, 2γ-2im, -2γ-2im, -γ-1.5im, -0.5-2im, -2γ+2im, -2γ, 0.1+0.2im, 0.1+im, 0.1-im, 0.8+0.1im, -0.8-0.1im,
                        γ+0.1im, -γ+0.1im, γ-0.3im, -γ-0.3im, γ+im, -γ+im, -γ-im, γ-im, γ-0.0im, -γ-0.0im)
                L = logkernelrectangle(z, n, γ)
                for j = 0:n-1, k=0:n-j-1
                    q = quadgk(s -> quadgk(t -> iszero(z-(s+im*t)) ? 0.0+0.0im : log(z-(s+im*t)) * legendrep(k,s/γ) * legendrep(j,t), -1, 1, atol=1E-12)[1], -γ, γ, atol=1E-12)[1]
                    @test q ≈ L[k+1,j+1] atol=1E-10
                end
            end
        end
    end

    @testset "BigFloat" begin
        setprecision(2000) do
            for γ in big.((1.3, 0.6))
                P = Legendre{BigFloat}()
                n = 150
                M = Diagonal((P'P)[1:n,1:n])
                for z in (2+γ, 2im*one(γ), 2+γ+2im, 2+γ-2im, -2-γ-2im)
                    @test P[big(0.), 1:n]' *  inv(M)  * logkernelrectangle(z, n, γ) * inv(M) * P[big(0.), 1:n] / γ ≈ log(z) atol=1E-40
                end
            end
            @test logkernelrectangle(Double64(0) + 0im, 10, Double64(1.3)) ≈ logkernelrectangle(big(0.0) + 0im, 10, big(1.3)) atol=1E-28
        end
    end
end

@testset "newtonian" begin
    @test newtonianrectangle(0.1+0.2im, 10) == newtonianrectangle([0.1,0.2], 10) == real(logkernelrectangle(0.1+0.2im, 10))
    @test newtonianrectangle(0.1+0.2im, 10, 1.3) == newtonianrectangle([0.1,0.2], 10, 1.3) == real(logkernelrectangle(0.1+0.2im, 10, 1.3))
end

# n = 4000
# z = 2.0
# @profview logkernelrectangle(z, n)

# x = range(-1,1,n); @time log.(abs.(x .+ im .* x'));