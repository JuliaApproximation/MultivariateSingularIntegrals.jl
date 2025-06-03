# MultivariateSingularIntegrals.jl
 Julia package for computing multivariate singular integrals


The function `newtoniansquare([x,y], n)` computes a matrix of  Newtonian potentials of Legendre polynomials on the unit square $Ω := [-1,1]^2$. That is it computes
```math
\begin{pmatrix}
L_{00}(𝐱) & \cdots & L_{0,p-1}(𝐱) & L_{0p}(𝐱) \\
L_{10}(𝐱) & \cdots & L_{1,p-1}(𝐱) \\
\vdots & \iddots \\
L_{p0}(𝐱)
\end{pmatrix}
```
where $𝐱 = [x,y]$ is any point in $ℝ²$ (including on or near the unit square $Ω$) and 
```math
L_{k,j}(𝐱) := ∬_Ω P_k(s) P_j(t) \log(\| [s,t] - [x,y] \|) \, ds \, dt
```
where $P_k$ and $P_j$ are the Legendre polynomials of degree $k$ and $j$, respectively.
