using BenchmarkTools
using Test
using GLMakie
using LaTeXStrings

include("naiv_inter.jl")
include("lagrange_inter.jl")
include("newton_inter.jl")

function test_all()
    
    # naiv_inter
    println("testing naiv_inter()")
    println("test a: ", @test naiv_inter([0.0, 1, 2, 3], [-5.0, -6, -1, 16], [0.5, 1.5, 2.5]) ≈ [−5.8750, −4.6250, 5.6250])
    println("test b: ", @test naiv_inter([0, 0.5, 1, 1.5, 2, 2.5, 3], [-5.0, -6, -1, 16, 10, 20, 9], [3.5, 5]) ≈ [−551.0, −34305])
    
    # lagrange_inter
    println("testing lagrange_inter()")
    println("test a: ", @test lagrange_inter([0.0, 1, 2, 3], [-5.0, -6, -1, 16], [0.5, 1.5, 2.5]) ≈ [−5.8750, −4.6250, 5.6250])
    println("test b: ", @test lagrange_inter([0, 0.5, 1, 1.5, 2, 2.5, 3], [-5.0, -6, -1, 16, 10, 20, 9], [3.5, 5]) ≈ [−551.0, −34305])

    # newton_inter
    println("testing newton_inter()")
    println("test a: ", @test newton_inter([0.0, 1, 2, 3], [-5.0, -6, -1, 16], [0.5, 1.5, 2.5]) ≈ [−5.8750, −4.6250, 5.6250])
    println("test b: ", @test newton_inter([0, 0.5, 1, 1.5, 2, 2.5, 3], [-5.0, -6, -1, 16, 10, 20, 9], [3.5, 5]) ≈ [−551.0, −34305])

    return nothing
end

function plot_interp_all(n, pos_in_figure, func, fname)

    x = range(0, 2, length=n)
    y = func.(x)
    
    m = 2000 # "sehr hohe Auflösung"
    x_high_res = range(-10, 12, length=m)

    axis = Axis(pos_in_figure,
                title = "N = $n",
                ylabel = "y",
                xlabel = "x",
                limits = ((nothing, nothing), (-5, 5))) # specify y axis-range (keep x axis-range untouched)

    # func(xi)
    lines!(axis, x_high_res, func.(x_high_res),
           label=latexstring("\$ $(fname)(x)\$"),
           linestyle=:dash,
           linewidth=2)

    # Naiv
    lines!(axis, x_high_res, naiv_inter(x, y, x_high_res),
           label=L"p_{naiv}(x)",
           linewidth=3)

    # Lagrange
    lines!(axis, x_high_res, lagrange_inter(x, y, x_high_res),
           label=L"p_{lagrange}(x)",
           linewidth=2)

    # Newton
    lines!(axis, x_high_res, newton_inter(x, y, x_high_res),
           label=L"p_{newton}(x)",
           linewidth=3)

    
    scatter!(axis, x, y, label="Stützstellen")

    axislegend(position = :rb, labelsize=22)

    return nothing
end

function plot_interp_all(func; fname="f")

    # create a global figure (we plot stuff into the figure)
    figure = Figure(fontsize=22)

    plot_interp_all(5, figure[1,1], func, fname)
    plot_interp_all(15, figure[1,2], func, fname)
    plot_interp_all(25, figure[2,1], func, fname)
    plot_interp_all(50, figure[2,2], func, fname)

    display(GLMakie.Screen(), figure)

    return nothing
end

function plot_interp_error_all(n, pos_in_figure, func, fname)

    x = range(0, 2, length=n)
    y = func.(x)
    
    m = 2000 # "sehr hohe Auflösung"
    x_high_res = range(-5, 7, length=m)

    axis = Axis(pos_in_figure,
                title = latexstring("\$\\textrm{Fehler} \\quad \\left| p(x) - $(fname)(x) \\right| \\quad \\textrm{bei N = $n Stützstellen}\$"),
                ylabel = latexstring("\$\\textrm{Fehler} \\quad \\left| p(x) - $(fname)(x) \\right|\$"),
                xlabel = "x",
                yscale = log10)

    naiv_error = abs.(naiv_inter(x, y, x_high_res) - func.(x_high_res))
    
    lines!(axis, x_high_res, naiv_error,
           label="Naiv",
           linewidth=3,
           color=naiv_error,
           colormap=:rainbow1,
           colorrange=(0, 1e-4))

    lagrange_error = abs.(lagrange_inter(x, y, x_high_res) - func.(x_high_res))

    lines!(axis, x_high_res, lagrange_error,
           label="Lagrange",
           linewidth=3,
           color=lagrange_error,
           colormap=:rainbow1,
           linestyle=(:dot, :dense),
           colorrange=(0, 1e-4))

    newton_error = abs.(newton_inter(x, y, x_high_res) - func.(x_high_res))

    lines!(axis, x_high_res, newton_error,
           label="Newton",
           linewidth=3,
           color=newton_error,
           colormap=:rainbow1,
           linestyle=(:dash, :dense),
           colorrange=(0, 1e-4))
    
    vlines!(x, label="Stützstellen", color=:grey)

    axislegend(position = :rb, size = 2)

    return nothing
end

function plot_interp_error_all(func; fname="f")

    figure = Figure(fontsize=22)

    plot_interp_error_all(5, figure[1,1], func, fname)
    plot_interp_error_all(15, figure[1,2], func, fname)
    plot_interp_error_all(25, figure[2,1], func, fname)
    plot_interp_error_all(50, figure[2,2], func, fname)

    display(GLMakie.Screen(), figure)

    return nothing
end

function plot_interp_error_scaling_all(func; fname="f")

    figure = Figure(fontsize=22)

    N = 25
    n_values::Vector{Int} = Int.(trunc.(1.2 .^ (4:N+3)))

    m = 2000 # "sehr hohe Auflösung"

    axis = Axis(figure[1,1],
                title = latexstring("\$\\textrm{Fehler} \\quad \\left| p(x) - $(fname)(x) \\right|_\\infty \\quad x \\in [0, 2] \\quad \\textrm{ beim Erhöhen der Anzahl an Stützstellen N und Auflösung M = $m}\$"),
                ylabel = latexstring("\$\\textrm{Fehler} \\quad \\left| p(x) - $(fname)(x) \\right|_\\infty \\quad x \\in [0, 2]\$"),
                xlabel = L"\textrm{Stützstellen N}",
                xscale = log10,
                yscale = log10,
                xticks = n_values)

    naiv_max_error = Vector{Float64}(undef, N)
    lagrange_max_error = Vector{Float64}(undef, N)
    newton_max_error = Vector{Float64}(undef, N)

    for i in 1:N

        x = range(0, 2, length=n_values[i])
        y = func.(x)
        
        x_high_res = range(0, 2, length=m)

        naiv_max_error[i] = maximum(abs.(naiv_inter(x, y, x_high_res) - func.(x_high_res)))
        lagrange_max_error[i] = maximum(abs.(lagrange_inter(x, y, x_high_res) - func.(x_high_res)))
        newton_max_error[i] = maximum(abs.(newton_inter(x, y, x_high_res) - func.(x_high_res)))
    end
    
    lines!(axis, n_values, naiv_max_error,
           label="Naiv",
           linewidth=3,
           color=naiv_max_error,
           colormap=:rainbow1,
           colorrange=(0, 1e-4))

    lines!(axis, n_values, lagrange_max_error,
           label="Lagrange",
           linewidth=3,
           color=lagrange_max_error,
           colormap=:rainbow1,
           linestyle=(:dot, :dense),
           colorrange=(0, 1e-4))

    lines!(axis, n_values, newton_max_error,
           label="Newton",
           linewidth=3,
           color=newton_max_error,
           colormap=:rainbow1,
           linestyle=(:dash, :dense),
           colorrange=(0, 1e-4))

    axislegend(position = :rb)

    display(GLMakie.Screen(), figure)
    
    return nothing
end

function plot_benchmark_all(m, pos_in_figure, func)
    
    naiv_inter_time_ms = zeros(Float64, 4)
    lagrange_inter_time_ms = zeros(Float64, 4)
    newton_inter_time_ms = zeros(Float64, 4)

    resolutions = [5, 15, 25, 50]

    for i in 1:4
        x = range(0, 2, length=resolutions[i])
        y = func.(x)
        
        x_high_res = range(-5, 7, length=m)

        bench_naiv_inter = @benchmark naiv_inter($x, $y, $x_high_res)
        naiv_inter_time_ms[i] = median(bench_naiv_inter).time * 1e-3

        bench_lagrange_inter = @benchmark lagrange_inter($x, $y, $x_high_res)
        lagrange_inter_time_ms[i] = median(bench_lagrange_inter).time * 1e-3

        bench_newton_inter = @benchmark newton_inter($x, $y, $x_high_res)
        newton_inter_time_ms[i] = median(bench_newton_inter).time * 1e-3
    end

    colors = Makie.wong_colors()

    axis = Axis(pos_in_figure,
                title = "Auswertung des Polynoms an M = $m Stellen",
                ylabel = "μs",
                xticks = (1:3, ["Naiv", "Lagrange", "Newton"]))
    
    barplot!(axis,
             [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
             [naiv_inter_time_ms..., lagrange_inter_time_ms..., newton_inter_time_ms...],
             dodge = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4],
             bar_labels = :y,
             label_formatter = x-> "$(round(x, sigdigits=3)) μs",
             color = colors[[1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]])

    return nothing
end

function plot_benchmark_all(func)
    
    figure = Figure(fontsize=22)

    plot_benchmark_all(2, figure[1,1], func)
    plot_benchmark_all(2000, figure[1,2], func)

    colors = Makie.wong_colors()

    Legend(figure[1,3],
           [PolyElement(polycolor = colors[i]) for i in 1:4],
           ["N = 5", "N = 15", "N = 25", "N = 50"])

    display(GLMakie.Screen(), figure)

    return nothing
end
