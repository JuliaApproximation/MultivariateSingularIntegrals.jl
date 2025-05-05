"""
    stieltjessquare(z, n)

computes the matrix with entries ``∫_{-1}^1∫_{-1}^1 P_k(s)P_j(t)/(z-(s+im*t)) dsdt`` up to total degree ``n``.
The bottom right of the returned matrix is zero. For a square truncation compute `stieltjessquare(z,2n-1)[1:n,1:n]`.
"""

function stieltjessquare(z, n)
    T = complex(float(eltype(z)))
    stieltjessquare!(zeros(T,n,n), z)
end

function stieltjessquare_populatefirstcolumn!(A, z)
    M_p = imlogkernel_vec(n, z+1)
    M_m = imlogkernel_vec(n, z-1)
    A[:,1] .= M_p .- M_m
end

function stieltjessquare_populatefirstrow!(A, z)
    M_p = imlogkernel_vec(n, 1-im*z)
    M_m = imlogkernel_vec(n, -1-im*z)
    m = size(A,1)
    A[1,:] .= (-1) .^ (1:m) .* im .* (M_p .- M_m)
end


function stieltjessquare!(A::AbstractMatrix{T}, z) where T
    m,n = size(A)
    @assert m  == n

    stieltjessquare_populatefirstcolumn!(A, z, F_1, F_2)
    stieltjessquare_populatefirstrow!(A, z, F_1, F_2)

    # 2nd row/column
    for k = 1:m-2
        A[k+1,2] = im* (k*A[k,1]/(2k+1)  + (k+1)*A[k+2,1]/(2k+1) - z*A[k+1,1])
    end

    for j = 2:n-2
        A[2,j+1] = z*A[1,j+1] - im*(j*A[1,j]/(2j+1) + (j+1)*A[1,k+2]/(2j+1))
    end
    # remaining
    for ℓ = 1:((n-1)÷2-1)
        for k = ℓ+1:n-(ℓ+2)
            A[k+1,ℓ+2] = (im*(2ℓ+1)*(k*A[k,ℓ+1]/(2k+1)  + (k+1)*A[k+2,ℓ+1]/(2k+1) - z*A[k+1,ℓ+1]) - j*A[k+1,ℓ+1])/(ℓ+1)
        end
        for j = ℓ+2:n-(ℓ+2)
            A[ℓ+2,j+1] = (z*(2k+1)*A[ℓ+1,j+1] - im*(2k+1)*(j*A[ℓ+1,j]/(2j+1) + (j+1)*A[ℓ+1,k+2]/(2j+1)) - k*A[k,ℓ+1])/(k+1)
        end
    end
    A
end