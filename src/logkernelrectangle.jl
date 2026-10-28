"""
    zlog(z)

implements ``z*log(z)``.
"""
zlog(z) = iszero(z) ? zero(z) : z*log(z)

"""
    zlogm(z)

implements ``z*log(z)`` but taking the other choice on the branch cut.
"""
zlogm(z) = iszero(imag(z)) && real(z) < 0 ? -zlog(abs(z)) - im*π*z : zlog(z)

L0(z) = zlog(1 + complex(z)) - zlog(complex(z)-1) - 2one(z)
L1(z, r0=L0(z)) = (z+1)*r0/2 + 1 - zlog(complex(z)+1)

function m_const(k, z)
    (x,y) = reim(z)
    T = float(real(typeof(z)))
    im*convert(T,π) * if k == 0
        if x > 0 || y ≥ 1 || (x == 0 && !signbit(x) && signbit(y)) # need to treat im*y+0.0 differently
            convert(T, 1)
        elseif y ≤ -1 # x ≤ 0
            convert(T, -3)
        else #  x ≤ 0 and -1 < y < 1
            ((1 + y) - 3*(1-y))/2
        end
    elseif -1 ≤ y ≤ 1 && x ≤ 0 && !(x == 0 && !signbit(x) && signbit(y)) # k ≠ 0
        -Weighted(Jacobi{T}(1,1))[y,k]/k
    else
        zero(T)
    end
end

M0(z) = L0(-im*float(z)) + m_const(0, float(z))
M1(z, r0=L0(-im*float(z))) = L1(-im*float(z), r0) + m_const(1, float(z))


"""
    L₀₀(z, γ=1) = ∫_{-γ}^γ∫_{-1}^1 log(z-(s+im*t))dtds
"""
L₀₀(z, γ=1) = (γ-z)*M0(z-γ) + im*M1(z-γ) + (γ+z)*M0(z+γ) - im*M1(z+γ) - 4γ

function m_const_vec(n, z)
    (x,y) = reim(z)
    T = float(real(typeof(z)))
    if x > 0 || y ≥ 1 || (x == 0 && !signbit(x) && signbit(y)) # need to treat im*y+0.0 differently
        [im*convert(T,π); Zeros{T}(n-1)]
    elseif y ≤ -1
        [-3im*convert(T,π); Zeros{T}(n-1)]
    else # x ≤ 0 && -1 ≤ y ≤ 1
        r = Weighted(Jacobi{T}(1,1))[y,1:n-1]
        [((1 + y) - 3*(1-y))/2*im*convert(T,π); -im .* convert(T,π) .* r ./ (1:(n-1))]
    end
end



function rec_rhs_1!(F::AbstractMatrix{T}, z, γ=one(real(T))) where T
    m,n = size(F)
    x,y = reim(z)
    πT = convert(T, π)
    fill!(F, zero(T))
    if x < -γ && -1 < y < 1
        C_y = Ultraspherical{T}(-1/2)[y,2:n+1]
        F[1,:] .=  (-4im*πT * x) .* C_y
    elseif x < γ && -1 < y < 1
        C_y = Ultraspherical{T}(-1/2)[y,2:n+1]
        F[1,:] .=  (-2im*πT * (x-γ)) .* C_y
    end

    F[1,1] += zlog(z-γ-im) + zlogm(z-γ+im) + zlog(z+γ-im) + zlogm(z+γ+im)
    F[1,2] -= 4im/convert(T,3)

    M₋ = imlogkernel_vec(n+1, z-γ)
    M₊ = imlogkernel_vec(n+1, z+γ)

    F[1,1] += im*(M₊[2] + M₋[2])
    for j = 1:n-1
        F[1,j+1] += im*(M₊[j+2] + M₋[j+2] - M₊[j] - M₋[j])/(2j+1)
    end
    F[1,:] .*= γ
    F[2,1] -= 4γ^2/convert(T,3)
    F
end



function rec_rhs_2!(F::AbstractMatrix{T}, z, γ=one(real(T))) where T
    m,n = size(F)
    x,y = reim(z)
    πT = convert(T, π)
    if -γ < x < γ && -1 ≤ y < 1
        C_x = Ultraspherical{T}(-3/2)[x/γ,3:m+2]
        C_y = Ultraspherical{T}(-1/2)[y,2:n+1]
        F .= (2im*πT*γ^2) .* (C_x .* C_y') ./ 3
        F[1,:] .-= (2im*πT*γ) .* x .* C_y
        F[2,:] .+= (2im*πT*γ^2/3) .* C_y
    elseif x ≤ -γ && -1 ≤ y < 1
        fill!(F, zero(T))
        C_y = Ultraspherical{T}(-1/2)[y,2:n+1]
        F[1,:] .= (-4im*πT*γ) .* x .* C_y
        F[2,:] .= (4im*πT*γ^2/3) .*  C_y
    else
        fill!(F, zero(T))
    end

    # integrals over s ∈ [-γ,γ] are rescaled to [-1,1], with z ± im ↦ (z ± im)/γ
    z₋, z₊ = (z-im)/γ, (z+im)/γ
    L₋ = complexlogkernel(Legendre{T}(), z₋)[1:m+1]
    L₊ = complexlogkernel(Legendre{T}(), z₊)[1:m+1]

    F[1,1] += γ^2*(L₋[2] + L₊[2] + zlog(z₋-1) + zlog(z₋+1) + zlog(z₊-1) + zlog(z₊+1)) + 4γ*z*log(γ)
    F[1,2] -= 4im*γ/convert(T,3)
    for k = 1:m-1
        F[k+1,1] += γ^2*(L₊[k+2] + L₋[k+2] - L₊[k] - L₋[k])/(2k+1)
    end
    F[2,1] -= 4γ^2*(1+log(γ))/convert(T,3)
    F
end


function logkernelrectangle_populatefirstcolumn!(A, z, F_1, F_2, γ=1)
    A[2,1] = z * A[1,1]/(3γ) + (F_2[1,1] - 2F_1[1,1])/(3γ)
    for k = 1:size(A,1)-2
        A[k+2,1] = (2k+1) * z * A[k+1,1]/((k+3)γ) - (k-2)*A[k,1]/(k+3) + (2k+1) * (F_2[k+1,1] - 2F_1[k+1,1])/((k+3)γ)
    end
    A
end

function logkernelrectangle_populatefirstrow!(A, z, F_1, F_2)
    A[1,2] = z*A[1,1]/(3im) + (F_1[1,1]-2F_2[1,1])/(3im)
    for j = 1:size(A,2)-2
        A[1,j+2] = -im * (2j+1) * z * A[1,j+1]/(j+3) - (j-2)*A[1,j]/(j+3) - im * (2j+1) * (F_1[1,j+1] - 2F_2[1,j+1])/(j+3)
    end
    A
end


"""
    logkernelrectangle(z, n, γ=1)

computes the matrix with entries ``∫_{-γ}^γ∫_{-1}^1 log(z-(s+im*t))P_k(s/γ)P_j(t)dtds`` up to total degree ``n``.
The default `γ = 1` corresponds to the unit square.
The bottom right of the returned matrix is zero. For a square truncation compute `logkernelrectangle(z,2n-1,γ)[1:n,1:n]`.
"""
function logkernelrectangle(z, n, γ=1)
    T = complex(float(promote_type(eltype(z), typeof(γ))))
    logkernelrectangle!(zeros(T,n,n), z, zeros(T,n,n), zeros(T,n,n), γ)
end


function logkernelrectangle!(A::AbstractMatrix{T}, z, F_1, F_2, γ=one(real(T))) where T
    m,n = size(A)
    @assert m  == n
    γ = convert(real(T), γ)
    @assert γ > 0
    n == 0 && return A
    A[1,1] = L₀₀(z, γ)
    n == 1 && return A
    rec_rhs_1!(F_1, z, γ)
    rec_rhs_2!(F_2, z, γ)
    logkernelrectangle_populatefirstcolumn!(A, z, F_1, F_2, γ)
    logkernelrectangle_populatefirstrow!(A, z, F_1, F_2)

    F = F_1 # reuse the memory
    F .= F_2 .- F_1

    # 2nd row/column

    for k = 1:m-2
        A[k+1,2] = im*(F[k+1,1] + γ*(A[k,1] - A[k+2,1])/(2k+1))
    end

    for j = 2:n-2
        A[2,j+1] = (F[1,j+1] + im*(A[1,j+2] - A[1,j])/(2j+1))/γ
    end

    @inbounds for ℓ = 1:((n-1)÷2-1)
        for k = ℓ+1:n-(ℓ+2)
            A[k+1,ℓ+2] = (2ℓ+1)*im*(F[k+1,ℓ+1] + γ*(A[k,ℓ+1] - A[k+2,ℓ+1])/(2k+1)) + A[k+1,ℓ]
        end
        for j = ℓ+2:n-(ℓ+2)
            A[ℓ+2,j+1] = (2ℓ+1) * (F[ℓ+1,j+1] + im*(A[ℓ+1,j+2] - A[ℓ+1,j])/(2j+1))/γ + A[ℓ,j+1]
        end
    end
    A
end
