using Pkg
using Plots
using StatsBase
using Distributions
using StatsPlots
using LinearAlgebra
using Random
using Plots.PlotMeasures

Pkg.activate(".")
using CanEvolve

function plotter(sim_object; bins::Int64 = 100)
    VAF = sim_object[2][1]
    scaler = fit(Histogram, VAF, nbins = bins)
    scaler = normalize(scaler, mode = :pdf)
    plot(scaler, fillcolor = "#008EA0", linecolor = "#008EA0", linealpha = 0.5, legend = false, xlims = (0,1), ytickfont = font(8, "Helvetica"), xtickfont = font(8, "Helvetica"), xtickfontsize=8,ytickfontsize=8)
    if sum(sim_object[3][1]) != 0
        vline!(sim_object[3][1], color = "#FF6F00", alpha = 1, linewidth = 3)
    end
    xlabel!("Variant allele frequency (VAF)")
    ylabel!("Density")
    #density!(VAF, linecolor = RGBA(249/255,101/255,103/255,0), linewidth = 5, linealpha = 0.95)
end

anim = Animation()
function plotter_two(x, y)
    gr(display_type=:inline)
    pl = plot(plotVAF2(x, bins = 100), plotVAF2(y, bins = 100), size=(1000,400), bottom_margin = 10mm)
end
for j in 1:20
    p, n = autoSimulation_fixed(log(2), 0.0, u = 100, depth = 120, rho = 0.001, lower_cutoff = 0.2, upper_cutoff = 0.35)
    push!(p, j, plotter_two(p, n))
    frame(anim)
end

gif(anim, "img/autosimulation.gif", fps=1)
