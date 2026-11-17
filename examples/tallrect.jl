using MultivariateSingularIntegrals

P = Legendre()
a,b,c,d = 0,1,0,3; # (a,b) × (c,d)

f = (x,y) -> cos(x * exp(y));

𝐱 = [0.1, 0.2]; # point inside the rectangle
@time qgk_approx = quadgk(s -> quadgk(t -> f(s, t) * log(norm(𝐱 - [s,t])), c, 𝐱[2], d)[1], a, 𝐱[1], b)[1]

h, γ, 𝐱₀ = (d-c)/2, (b-a)/(d-c), [(a+b)/2, (c+d)/2];
p = 60;
x = ChebyshevGrid{1}(p)
y = x'
F = f.((a+b)/2 .+ (b-a)/2 .* x, (c+d)/2 .+ (d-c)/2 .* y); # evaluate f on the mapped grid
C = plan_transform(P, (p, p)) * F; # Legendre coefficients c_kj




@test sum(C[k,j]*P[0.1,k]*P[0.2,j] for k=1:p, j=1:p) ≈ f((a+b)/2 .+ (b-a)/2 .* 0.1, (c+d)/2 .+ (d-c)/2 .* 0.2)

@time N = newtonianrectangle((𝐱 - 𝐱₀)/h, p, γ);

@time N = setprecision(5000) do
    Float64.(newtonianrectangle(big.((𝐱 - 𝐱₀)/h), p, big(γ)))
end

@test h^2 * (dot(N, C) + 4γ*log(h)*C[1,1]) ≈ qgk_approx


γ = 0.5 # [-γ,γ] × [-1,1]

𝐱 = [1.1, 0.2]
@test quadgk(s -> quadgk(t -> log(norm(𝐱 - [s,t])), -1, 1)[1], -γ, γ)[1] ≈ newtonianrectangle(𝐱, p, γ)[1]

𝐱 = [0.1, 0.2]
p = 30; norm(newtonianrectangle(𝐱, p, γ) - newtonianrectangle(big.(𝐱), p, big(γ)))
using Chairmarks
[1:20,1]


γ = 0.01
𝐱 = [0.001, 0.2]
p = 10; norm(newtonianrectangle(𝐱, p, γ) - newtonianrectangle(big.(𝐱), p, big(γ)))

newtonianrectangle(Float64x2.(𝐱), p, Float64x2.(γ))


γ = 100
𝐱 = [0.001, 0.2]
p = 10; norm(newtonianrectangle(𝐱, p, γ) - newtonianrectangle(big.(𝐱), p, big(γ)))


