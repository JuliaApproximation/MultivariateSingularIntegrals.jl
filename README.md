# MultivariateSingularIntegrals.jl
 Julia package for computing multivariate singular integrals


The function `newtonianrectangle([x,y], p)` computes a matrix of  Newtonian potentials of Legendre polynomials on the unit square $Ω := [-1,1]^2$. That is it computes
```math
\begin{pmatrix}
L_{11}(𝐱) & ⋯ & L_{1,p-1}(𝐱) & L_{1p}(𝐱) \\
L_{21}(𝐱) & ⋯ & L_{2,p-1}(𝐱) \\
⋮ & ⋰ \\
L_{p1}(𝐱)
\end{pmatrix}
```
where $𝐱 = [x,y]$ is any point in $ℝ²$ (including on or near the unit square $Ω$) and 
```math
L_{k,j}(𝐱) := ∬_Ω P_{k-1}(s) P_{j-1}(t) \log(\| [s,t] - 𝐱 \|) \, ds \, dt
```
where $P_k$ and $P_j$ are the Legendre polynomials of degree $k$ and $j$, respectively.
More generally, `newtonianrectangle([x,y], p, γ)` computes the same quantities on the rectangle $Ω := [-γ,γ] × [-1,1]$
with $P_{k-1}(s)$ replaced by $P_{k-1}(s/γ)$.


Here we show an example of how to use the package for computing the Newtonian potential for $f(x,y) = \cos(x \exp(y))$,
which is faster than QuadGK.jl:
```julia
julia> using MultivariateSingularIntegrals, ClassicalOrthogonalPolynomials, LinearAlgebra, QuadGK

julia> f = (x,y) -> cos(x * exp(y));

julia> 𝐱 = [0.1, 1.1]; # point near the unit square

julia> @time qgk_approx = quadgk(s -> quadgk(t -> f(s, t) * log(norm(𝐱 - [s,t])), -1, 1)[1], -1, 1)[1]; # compute the integral using quadgk
  0.388834 seconds (804.74 k allocations: 35.045 MiB, 3.99% gc time, 98.86% compilation time)

julia> P = Legendre(); # Legendre polynomials

julia> p = 40; # degree of the Legendre polynomials

julia> x = ClassicalOrthogonalPolynomials.grid(P, p); # grid of points

julia> F = f.(x, x') # evaluate f on the grid;

julia> @time C = plan_transform(P, (p, p)) * F; # 2D Legendre coefficients
  0.002116 seconds (2.87 k allocations: 3.594 MiB)

julia> @time N = Float64.(newtonianrectangle(big.(𝐱), p)); # compute the Newtonian potential of Legendre polynomials, using BigFloat to avoid numerical issues
  0.003480 seconds (80.94 k allocations: 4.119 MiB)

julia> @test dot(N, C) ≈ qgk_approx
Test Passed
```
The package continues to work inside the unit square:
```julia
julia> 𝐱 = [0.1, 0.2];

julia> @time N = Float64.(newtonianrectangle(big.(𝐱), p));
  0.005985 seconds (165.52 k allocations: 8.324 MiB)

julia> dot(N, C)
-1.2555132824835833
```

The package also supports complex logarithmic integrals with kernel $\log(z-(s+{\rm i}t))$ via `logkernelrectangle`
and Stieltjes integrals with kernel $1/(z-(s+{\rm i}t))$ via `stieltjessquare`. The real and imaginary parts of Stieltjes integrals containing the gradient of the Newtonian potential.
