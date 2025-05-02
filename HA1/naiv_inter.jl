
function naiv_inter(x, y, xi)
    @assert(length(x) == length(y))
    size = length(x)

    #          i
    #       ⎡     ⎤
    # A  =  ⎢     ⎥ j
    #       ⎣     ⎦
    
    A = Matrix{Float64}(undef, size, size)

    # Füllen der Matrix mit m_i(x_j) wobei m_i die i-te monomiale Basis ist.
    for i in 1:size
        for j in 1:size
            A[j, i] = x[j] ^ (i - 1)
        end
    end

    # Lösen des LGS
    coeffs = inv(A) * y

    result = Vector{Float64}(undef, length(xi))
    for i in 1:length(xi)
        result[i] = 0
        for j in 1:size
            result[i] += (xi[i] ^ (j - 1)) * coeffs[j]
        end
    end
    
    return result
end
