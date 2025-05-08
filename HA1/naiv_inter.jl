using LinearSolve

function naiv_inter(x, y, xi)
    @assert(length(x) == length(y))
    N = length(x) # Anzahl Stützstellen
    M = length(xi) # Anzahl Auswertungspunkte

    #          i
    #       ⎡     ⎤
    # A  =  ⎢     ⎥ j
    #       ⎣     ⎦
    
    A = Matrix{Float64}(undef, N, N)

    # Füllen der Matrix mit m_i(x_j) wobei m_i die i-te monomiale Basis ist.
    for i in 1:N
        for j in 1:N
            A[j, i] = x[j] ^ (i - 1)
        end
    end

    # Lösen des LGS
    # coeffs = inv(A) * y
    coeffs = A\y

    yi = zeros(Float64, M)
    for i in 1:M
        for j in 1:N
            yi[i] += (xi[i] ^ (j - 1)) * coeffs[j]
        end
    end
    
    return yi
end
