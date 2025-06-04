using LinearAlgebra
using Test
using GLMakie
using LaTeXStrings

function gauss_quadrature_n(f, a, b, n)::Float64

    # Berechnen der Stützstellen
    A = zeros(n, n)
    for k in 1:(n - 1)
        β = k / √(4 * k^2 - 1);
        A[k, k + 1] = β
        A[k + 1, k] = β
    end

    xs, Z = eigen(SymTridiagonal(A)) # Small speedup by specifying matrix-form 'SymTridiagonal'
    
    # Berechnen der Gewichte
    ws = 2 * Z[1, :] .^ 2

    # Integrieren
    sum::Float64 = 0
    for i in 1:n
        sum += ws[i] * f(((b - a) / 2) * xs[i] + (a + b) / 2)
    end

    return ((b - a) / 2) * sum
end

function gauss_quadrature_tol_naiv(f, a, b, tol)

    n::Int32 = 0
    n_vec = Vector{Int32}()
    
    Q_prev = gauss_quadrature_n(f, a, b, n+=1)
    Q = gauss_quadrature_n(f, a, b, n+=1)
    push!(n_vec, n)

    E = abs(Q - Q_prev)

    invalid_E = isnan(E) || isinf(E)
    
    while E > tol || invalid_E

        n += 1

        # updating Q
        Q_prev = Q
        Q = gauss_quadrature_n(f, a, b, n)

        # updating E
        E = abs(Q - Q_prev)
        
        invalid_E = isnan(E) || isinf(E)

        push!(n_vec, n)
    end

    return Q, n
end

function gauss_quadrature_tol(f, a, b, tol)

    n::Int32 = 0
    Δn = 1
    n_vec = Vector{Int32}()

    # Erste 2 Iterationen, um E berechnen zu können
    Q_vec = Vector{Float64}()
    Q_prev = gauss_quadrature_n(f, a, b, n+=1)
    push!(Q_vec, Q_prev)
    push!(n_vec, n)
    Q = gauss_quadrature_n(f, a, b, n+=1)
    push!(Q_vec, Q)
    push!(n_vec, n)

    E = abs(Q - Q_prev) / Δn

    lnΔE = -20 # So gewählt, dass Δn am Anfang klein ist.
    lnΔE_mov_avg = lnΔE
    E_vec = Vector{Float64}()
    push!(E_vec, NaN) # Make this the same size as 'n_vec'
    push!(E_vec, E)

    println("n: ", n, " Δn: ", Δn, " E: ", E , " lnΔE: ", lnΔE, " lnΔE_mov_avg: ", lnΔE_mov_avg)

    invalid_E = isnan(E) || isinf(E)
    
    while E > tol || invalid_E

        # updating n
        prev_n = n
        if invalid_E
            n += 1
        else
            lnΔE_mov_avg = lnΔE_mov_avg * 0.9 + lnΔE * 0.1
            n = ceil(log(tol / E) * Δn / -abs(lnΔE_mov_avg)) + n
        end
        Δn = n - prev_n
        push!(n_vec, n)

        # updating Q
        Q_prev = Q
        Q = gauss_quadrature_n(f, a, b, n)
        push!(Q_vec, Q)

        # updating E
        prev_E = E
        E = abs(Q - Q_prev) / Δn
        push!(E_vec, E)

        println("n: ", n, " Δn: ", Δn, " E: ", E, " lnΔE: ", lnΔE, " lnΔE_mov_avg: ", lnΔE_mov_avg)
        
        lnΔE = log(E / prev_E)
        invalid_E = isnan(lnΔE) || isinf(lnΔE)
    end

    return Q, n, n_vec, Q_vec, E_vec
end

function test()

    figure = Figure(fontsize=22)
    
    axis = Axis(figure[1,1],
                title = latexstring("\$\\textrm{Fehler} \\quad \\left| Q_n - ∫ f(x)dx \\right| \\quad \\quad \\textrm{ In Abhängigkeit von n, für verschiedene Funktionen}\$"),
                ylabel = latexstring("\$\\textrm{Fehler} \\quad \\left| Q_n - ∫ f(x)dx \\right| \\quad \$"),
                xlabel = L"\textrm{n}",
                yscale = log10)

    println("testing 'gauss_quadrature_n(f, a, b, n)':\n")
    
    println("[-1, 1] ∫ x^10 dx \t\t| num (n=6): ", gauss_quadrature_n((x -> x^10), -1, 1, 6), "\t| analytic: ", 2/11)
    println("[0, π] ∫ sin(x) dx \t\t| num (n=6): ", gauss_quadrature_n((x -> sin(x)), 0, pi, 6), "\t\t| analytic: ", 2)
    println("[-2, 3] ∫ 1 / (10^-2 + x^2) dx \t| num (n=1000): ", gauss_quadrature_n((x -> 1 / (10^-2 + x^2)), -2, 3, 1000), "\t| analytic: ", 10 * atan(20) + 10 * atan(30))

    println("\ntesting 'gauss_quadrature_tol(f, a, b, tol)':\n")

    error_tolerance = 1e-10

    println("[-1, 1] ∫ x^10 dx \t\t| target-error: ", error_tolerance)
    result, n, n_vec, Q_vec, E_vec = gauss_quadrature_tol((x -> x^10), -1, 1, error_tolerance)
    println("error: ", abs(result - 2/11), "\t | n: ", n, "\n")

    scatterlines!(axis, n_vec, abs.(Q_vec .- 2/11),
                  label="[-1, 1] ∫ x^10 dx | E=$(error_tolerance)",
                  linewidth=3, color=:red)
    scatterlines!(axis, n_vec, E_vec, linewidth=1, linestyle=:dash, color=:red)

    println("[0, π] ∫ sin(x) dx \t\t| target-error: ", error_tolerance)
    result, n, n_vec, Q_vec, E_vec = gauss_quadrature_tol((x -> sin(x)), 0, pi, error_tolerance)
    println("error: ", abs(result - 2), "\t | n: ", n, "\n")

    scatterlines!(axis, n_vec, abs.(Q_vec .- 2),
                  label="[0, π] ∫ sin(x) dx | E=$(error_tolerance)",
                  linewidth=3, color=:orange)
    scatterlines!(axis, n_vec, E_vec, linewidth=1, linestyle=:dash, color=:orange)

    println("[-2, 3] ∫ 1 / (10^-2 + x^2) dx \t| target-error: ", error_tolerance)
    result, n, n_vec, Q_vec, E_vec = gauss_quadrature_tol((x -> 1 / (10^-2 + x^2)), -2, 3, error_tolerance)
    println("error: ", abs(result - (10 * atan(20) + 10 * atan(30))), "\t | n: ", n, "\n")

    scatterlines!(axis, n_vec, abs.(Q_vec .- (10 * atan(20) + 10 * atan(30))),
                  label="[-2, 3] ∫ 1 / (10^-2 + x^2) dx | E=$(error_tolerance)",
                  linewidth=3, color=:blue)
    scatterlines!(axis, n_vec, E_vec, linewidth=1, linestyle=:dash, color=:blue)

    error_tolerance = 1e-7

    famous_cleo_integral(x) = (1/x) * sqrt((1+x) / (1-x)) * log((2x^2+2x+1) / (2x^2-2x+1))
    println("[-1, 1] ∫ (1/x) * √((1+x) / (1-x)) * log((2x^2+2x+1) / (2x^2-2x+1)) dx \t| target-error: ", error_tolerance)
    result, n, n_vec, Q_vec, E_vec = gauss_quadrature_tol(famous_cleo_integral, -1, 1, error_tolerance)
    println("error: ", abs(result - 4 * pi * acot(sqrt((1+√5)/2))), "\t | n: ", n, "\n")

    scatterlines!(axis, n_vec, abs.(Q_vec .- 4 * pi * acot(sqrt((1+√5)/2))),
                  label="[-1, 1] ∫ (1/x) * √((1+x) / (1-x)) * log((2x^2+2x+1) / (2x^2-2x+1)) dx | E=$(error_tolerance)",
                  linewidth=3, color=:green)
    scatterlines!(axis, n_vec, E_vec, linewidth=1, linestyle=:dash, color=:green)

    error_tolerance = 1e-5

    println("[0, 1] ∫ sin(1/x)/√x dx \t| target-error: ", error_tolerance)
    result, n, n_vec, Q_vec, E_vec = gauss_quadrature_tol((x -> sin(1/x) / √x), 0, 1, error_tolerance)
    println("error: ", abs(result - 0.5714732926457052), "\t | n: ", n, "\n")

    scatterlines!(axis, n_vec, abs.(Q_vec .- 0.5714732926457052),
                  label="[0, 1] ∫ sin(1/x)/√x dx | E=$(error_tolerance)",
                  linewidth=3, color=:purple)
    scatterlines!(axis, n_vec, E_vec, linewidth=1, linestyle=:dash, color=:purple)

    error_tolerance = 5e-4

    # Source: https://www-m3.ma.tum.de/bornemann/challengebook/Chapter1/chall_int.pdf
    println("[0, 1] ∫ 1/x * cos(1/x * ln(x)) dx \t| target-error: ", error_tolerance)
    result, n, n_vec, Q_vec, E_vec = gauss_quadrature_tol((x -> 1/x * cos(1/x * log(x))), 0, 1, error_tolerance)
    println("error: ", abs(result - 0.3233674316777787), "\t | n: ", n, "\n")

    scatterlines!(axis, n_vec, abs.(Q_vec .- 0.3233674316777787),
                  label="[0, 1] ∫ 1/x * cos(1/x * ln(x)) dx | E=$(error_tolerance)",
                  linewidth=3, color=:darkblue)
    scatterlines!(axis, n_vec, E_vec, linewidth=1, linestyle=:dash, color=:darkblue)

    
    axislegend(position = :rb)
    display(GLMakie.Screen(), figure)

    return nothing
end

function plot_required_n_for_tol()
    
    figure = Figure(fontsize=22)
    
    axis = Axis(figure[1,1],
                title = "Benötigtes n für gegebenen Fehler",
                ylabel = "Tolleranz Kriterium",
                xlabel = L"\textrm{n}",
                yscale = log10)


    error_tolerance = 1e-10

    n_vec1 = Vector{Int32}()
    n_vec2 = Vector{Int32}()
    n_vec3 = Vector{Int32}()

    tol_vec = Float64(10) .^ (-1:-1:-14)

    for error_tolerance in tol_vec
        _, n = gauss_quadrature_tol_naiv((x -> x^10), -1, 1, error_tolerance)
        push!(n_vec1, n)

        _, n = gauss_quadrature_tol_naiv((x -> sin(x)), 0, pi, error_tolerance)
        push!(n_vec2, n)

        _, n = gauss_quadrature_tol_naiv((x -> 1 / (10^-2 + x^2)), -2, 3, error_tolerance)
        push!(n_vec3, n)
    end

    scatterlines!(axis, n_vec1, tol_vec, label="[-1, 1] ∫ x^10 dx", linewidth=3, color=:red)
    scatterlines!(axis, n_vec2, tol_vec, label="[0, π] ∫ sin(x) dx", linewidth=3, color=:orange)
    scatterlines!(axis, n_vec3, tol_vec, label="[-2, 3] ∫ 1 / (10^-2 + x^2) dx", linewidth=3, color=:blue)
    
    axislegend(position = :rb)
    display(GLMakie.Screen(), figure)

    return nothing
end
