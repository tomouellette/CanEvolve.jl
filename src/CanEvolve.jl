module CanEvolve

using StatsBase
using Distributions
using LinearAlgebra
using DataFrames
using Random
using Plots
using StatsPlots
using PyCall

const np = PyNULL()
function __init__()
    copy!(np, pyimport_conda("numpy", "numpy"))
end

export

    # simulation.jl
    cell,
    mutate,
    primeSim,
    evolveTumour,
    BetaBin,
    sequence,
    detectableDrivers,
    processTumour,
    adjustProportion,
    mutationMetrics,
    simulation,
    full_synthetic_neutral,

    # automation.jl
    sampleParameters,
    autoSimulation,
    autoSimulation_fixed,
    autoSimulation_fixed_depth,
    plotVAF,
    frequencyThresholds,
    mapfrequency,
    uniformdensity,
    rescale,
    engineer,
    engineer_un

include("simulation.jl")
include("automation.jl")

end
