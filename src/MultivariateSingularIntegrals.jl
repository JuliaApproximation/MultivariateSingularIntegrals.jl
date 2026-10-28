module MultivariateSingularIntegrals
using LinearAlgebra, ClassicalOrthogonalPolynomials, SingularIntegrals, FillArrays
import Base: size

export newtonianrectangle, logkernelrectangle, logkernelrectangle!, stieltjessquare, stieltjessquare!


function imlogkernel_vec(n, z)
    T = float(real(typeof(z)))
    transpose(complexlogkernel(Legendre{T}(), -im*float(z)))[1:n] + m_const_vec(n, float(z))
end

"""
    newtonianrectangle(z, n, γ=1)

computes the matrix with entries ``∫_{-γ}^γ∫_{-1}^1 log(‖z-[s,t]‖)P_k(s/γ)P_j(t)dtds`` up to total degree ``n``,
where `z` is either a vector `[x,y]` or a complex number `x+im*y`. The default `γ = 1` corresponds to the unit square.
"""
newtonianrectangle(z::AbstractVector, n, γ=1) = real(logkernelrectangle(complex(z...), n, γ))
newtonianrectangle(z::Number, n, γ=1) = real(logkernelrectangle(z, n, γ))

include("stieltjessquare.jl")
include("logkernelrectangle.jl")


end # module MultivariateSingularIntegrals
