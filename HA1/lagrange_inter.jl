
function lagrange_inter(x, y, xi)
    @assert length(x)==length(y)
    
    yi = zeros(length(xi))
    
    #Setzen der Polynome über 3 Schleifen, aussetzen über if Statement
    for i in 1:length(xi)
        for j in 1:length(x)
            zähler = 1
            nenner = 1
            
            for k in 1:length(x)
                if k!=j
                    zähler *= (xi[i]-x[k])
                    nenner *= (x[j]-x[k])
                end
            end
        
            yi[i] += (zähler / nenner) * y[j]
        end
    end

    return yi
end
