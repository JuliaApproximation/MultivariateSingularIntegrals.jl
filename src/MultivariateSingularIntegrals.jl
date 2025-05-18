module MultivariateSingularIntegrals
using LinearAlgebra, ClassicalOrthogonalPolynomials, SingularIntegrals, FillArrays
import Base: size

export logkernelsquare, stieltjessquare, logkernelsquare!


function imlogkernel_vec(n, z)
    T = float(real(typeof(z)))
    transpose(complexlogkernel(Legendre{T}(), -im*float(z)))[1:n] + m_const_vec(n, float(z))
end

include("stieltjessquare.jl")
include("logkernelsquare.jl")


end # module MultivariateSingularIntegrals
