using BenchmarkTools
using Test
using GLMakie

include("naiv_inter.jl")
include("lagrange_inter.jl")
include("newton_inter.jl")

function test_all()
    
    # naiv_inter
    println("testing naiv_inter")
    println("test a: ", @test naiv_inter([0, 1, 2, 3], [-5, -6, -1, 16], [0.5, 1.5, 2.5]) ≈ [−5.8750, −4.6250, 5.6250])
    println("test b: ", @test naiv_inter([0, 0.5, 1, 1.5, 2, 2.5, 3], [-5, -6, -1, 16, 10, 20, 9], [3.5, 5]) ≈ [−551, −34305])
    
    # lagrange_inter
    println("testing lagrange_inter")
    println("test a: ", @test lagrange_inter([0, 1, 2, 3], [-5, -6, -1, 16], [0.5, 1.5, 2.5]) ≈ [−5.8750, −4.6250, 5.6250])
    println("test b: ", @test lagrange_inter([0, 0.5, 1, 1.5, 2, 2.5, 3], [-5, -6, -1, 16, 10, 20, 9], [3.5, 5]) ≈ [−551, −34305])

    # newton_inter
    println("testing newton_inter")
    println("test a: ", @test newton_inter([0, 1, 2, 3], [-5, -6, -1, 16], [0.5, 1.5, 2.5]) ≈ [−5.8750, −4.6250, 5.6250])
    println("test b: ", @test newton_inter([0, 0.5, 1, 1.5, 2, 2.5, 3], [-5, -6, -1, 16, 10, 20, 9], [3.5, 5]) ≈ [−551, −34305])

    return nothing
end

function interpolate_sin_and_plot_all(n, pos_in_figure)

    x = range(0, 2, length=n)
    y = sin.(x)
    
    m = 2000 # "sehr hohe Auflösung"
    x_high_res = range(-10, 12, length=m)

    axis = Axis(pos_in_figure,
                title = "n = $n",
                ylabel = "y",
                xlabel = "x",
                limits = ((nothing, nothing), (-5, 5))) # specify y axis-range (keep x axis-range untouched)

    # plot sin(xi)
    lines!(axis, x_high_res, sin.(x_high_res),
           label=L"sin(x)",
           linestyle=:dash,
           linewidth=2)

    lines!(axis, x_high_res, naiv_inter(x, y, x_high_res),
           label=L"p_{naiv}(x)",
           linewidth=2)

    lines!(axis, x_high_res, lagrange_inter(x, y, x_high_res),
           label=L"p_{lagrange}(x)",
           linewidth=2)

    lines!(axis, x_high_res, newton_inter(x, y, x_high_res),
           label=L"p_{newton}(x)",
           linewidth=2)

    scatter!(axis, x, y, label="Stützstellen")

    axislegend(position = :rb, labelsize=22)

    return nothing
end

function interpolate_sin_and_plot_all()

    # create a global figure (we plot stuff into the figure)
    figure = Figure(fontsize=22)

    interpolate_sin_and_plot_all(5, figure[1,1])
    interpolate_sin_and_plot_all(15, figure[1,2])
    interpolate_sin_and_plot_all(25, figure[2,1])
    interpolate_sin_and_plot_all(50, figure[2,2])

    return figure
end

function interpolate_sin_and_plot_error_all(n, pos_in_figure)

    x = range(0, 2, length=n)
    y = sin.(x)
    
    m = 2000 # "sehr hohe Auflösung"
    x_high_res = range(-5, 7, length=m)

    axis = Axis(pos_in_figure,
                title = "Fehler bei n = $n Stützstellen",
                ylabel = L"log_10(y)",
                xlabel = "x",
                yscale = log10)

    naiv_error = abs.(naiv_inter(x, y, x_high_res) - sin.(x_high_res))
    
    lines!(axis, x_high_res, naiv_error,
           label=L"\left| p_{naiv}(x) - sin(x) \right|",
           linewidth=3,
           color=naiv_error,
           colormap=:rainbow1,
           colorrange=(0, 1e-4))

    lagrange_error = abs.(lagrange_inter(x, y, x_high_res) - sin.(x_high_res))

    lines!(axis, x_high_res, lagrange_error,
           label=L"\left| p_{lagrange}(x) - sin(x) \right|",
           linewidth=3,
           color=lagrange_error,
           colormap=:rainbow1,
           linestyle=(:dot, :dense),
           colorrange=(0, 1e-4))

    newton_error = abs.(newton_inter(x, y, x_high_res) - sin.(x_high_res))

    lines!(axis, x_high_res, newton_error,
           label=L"\left| p_{newton}(x) - sin(x) \right|",
           linewidth=3,
           color=newton_error,
           colormap=:rainbow1,
           linestyle=(:dash, :dense),
           colorrange=(0, 1e-4))
    
    vlines!(x, label="Stützstellen", color=:grey)

    axislegend(position = :rb, size = 2)

    return nothing
end

function interpolate_sin_and_plot_error_all()

    figure = Figure(fontsize=22)

    interpolate_sin_and_plot_error_all(5, figure[1,1])
    interpolate_sin_and_plot_error_all(15, figure[1,2])
    interpolate_sin_and_plot_error_all(25, figure[2,1])
    interpolate_sin_and_plot_error_all(50, figure[2,2])

    return figure
end

function benchmark_and_plot_results_all(n, pos_in_figure)

    return nothing
end

function benchmark_and_plot_results_all()
    
    figure = Figure(fontsize=22)

    naiv_inter_time_ms = zeros(Float64, 4)
    lagrange_inter_time_ms = zeros(Float64, 4)
    newton_inter_time_ms = zeros(Float64, 4)

    resolutions = [5, 15, 25, 50]

    for i in 1:4
        x = range(0, 2, length=resolutions[i])
        y = sin.(x)
        
        m = 2000 # "sehr hohe Auflösung"
        x_high_res = range(-5, 7, length=m)

        bench_naiv_inter = @benchmark naiv_inter($x, $y, $x_high_res)
        naiv_inter_time_ms[i] = median(bench_naiv_inter).time * 1e-3

        bench_lagrange_inter = @benchmark lagrange_inter($x, $y, $x_high_res)
        lagrange_inter_time_ms[i] = median(bench_lagrange_inter).time * 1e-3

        bench_newton_inter = @benchmark newton_inter($x, $y, $x_high_res)
        newton_inter_time_ms[i] = median(bench_newton_inter).time * 1e-3
    end

    colors = Makie.wong_colors()

    axis = Axis(figure[1,1],
                title = "Performance Benchmark",
                ylabel = "μs",
                xticks = (1:3, ["Naiv", "Lagrange", "Newton"]))

    barplot!(axis,
             [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
             [naiv_inter_time_ms..., lagrange_inter_time_ms..., newton_inter_time_ms...],
             dodge = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4],
             bar_labels = :y,
             label_formatter = x-> "$(round(x, sigdigits=4)) μs",
             color = colors[[1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]])

    Legend(figure[1,2],
           [PolyElement(polycolor = colors[i]) for i in 1:4],
           ["N = 5", "N = 15", "N = 25", "N = 50"])

    return figure
end
