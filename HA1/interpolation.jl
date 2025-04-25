using BenchmarkTools
using Test

function naiv_inter(x, y, xi)
    @assert(length(x) == length(y))
    size = length(x)

    #      _  i _
    # A = |      | j
    #     |_    _|
    
    A = Matrix{Float64}(undef, size, size)

    # Füllen der Matrix mit m_i(x_j) wobei m_i die i-te monomiale Basis ist.
    for i in 1:size
        for j in 1:size
            A[j, i] = x[j] ^ (i - 1)
        end
    end

    A * a = y

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

function lagrange_inter(x, y, xi)

end

function newton_inter(x, y, xi)

end

function test_all()
    
    # naiv_inter
    println("testing naiv_inter")
    println("test a: ", @test naiv_inter([0, 1, 2, 3], [-5, -6, -1, 16], [0.5, 1.5, 2.5]) ≈ [−5.8750, −4.6250, 5.6250])
    println("test b: ", @test naiv_inter([0, 0.5, 1, 1.5, 2, 2.5, 3], [-5, -6, -1, 16, 10, 20, 9], [3.5, 5]) ≈ [−551, −34305])
    

end

function benchmark_all()

    # naiv_inter
    
end

println("hallo")
println("moin")
