
# function newton_inter_yunus(x, y, xi)
#     @assert(length(x) == length(y))
#     N = length(x) # Anzahl Stützstellen
#     M = length(xi) # Anzahl Auswertungspunkte

#     A = zeros(Float64, N, N)

#     div_diffs = zeros(Float64, N, N) # Dividierte Differenzen

#     for i in 1:N
#         div_diffs[i, 1] = y[i] # "Rekursionsanfang"
#         for j in 2:i
#             div_diffs[i, j] = (div_diffs[i, j-1] - div_diffs[i-1, j-1]) / (x[i] - x[i-j + 1])
#         end
#     end

#     return div_diffs

#     # Auswertung des Polynoms mit Horner-Schema:
#     # https://de.wikipedia.org/wiki/Polynominterpolation#Newtonscher_Algorithmus
#     yi = Vector{Float64}(undef, M)
#     for i in 1:M
#         b = div_diffs[N, N]
#         for j in (N-1):-1:1
#             b = b * (xi[i] - x[j]) + div_diffs[j, j]
#         end
#         yi[i] = b
#     end

#     return yi
# end

function newton_inter(x, y, xi)
    @assert length(x) == length(y)
    D = ones(length(x), length(x))
    D[:, 1] .= y
    a = ones(length(x))
    a[1] = D[1, 1]
    yi = zeros(length(xi))
    
    for k=2:length(x) 
       for j=2:length(x)
           if(j>=k) 
              D[j,k]=(D[j,k-1]-D[j-1,k-1])/(x[j]-x[j-k+1]); 
           end
       end 
       a[k]=D[k,k];
    end

    for i = 1:length(xi)
        temp = 1
        yi[i] = a[1]
        for n = 2:length(x)
            temp = temp * (xi[i] - x[n - 1])
            yi[i] = yi[i] + a[n] * temp
        end
    end

    return yi
end
