# ---------------------------------------------------------------------------
# automation.jl: Various functions for synthetic data generation and analysis
# ---------------------------------------------------------------------------

"""
Randomly sample parameters for tumour evolution simulation.

Input:
    - None

Returns:
    - integer: Mutation rate (per genome per doubling)
    - float: Probability of subclone/driver event (per division probability)
    - integer: Number of clonal mutations at t = 0
    - integer: Number of potential subclone/driver events during tumour growth
    - float: Mean of exponential selection coefficient distribution
    - integer: Mean sequencing depth
    - float: Sequencing overdispersion parameter for beta-binomial distribution
"""
function sampleParameters()
    u = rand(DiscreteUniform(1, 500), 1)[1]
    rho = rand(Uniform(0.0, 0.002), 1)[1]
    phi_b = rand(Uniform(0.001, 0.1), 1)[1] / u # Normalize by u so time of driver occurences are similar across different mutation rates
    nclonal = rand(DiscreteUniform(0, u*20), 1)[1]
    ndrivers = rand(DiscreteUniform(1, 3), 1)[1]
    lambda = rand(Uniform(1, 3), 1)[1]
    depth = rand(DiscreteUniform(50, 250), 1)[1]
    return u, phi_b, nclonal, ndrivers, lambda, depth, rho
end

"""
Generates synthetic sequenced tumour data for a population subject to positive selection and a population evolving neutrally.

Input:
    - b (float): Birth rate
    - d (float): Death rate
    - Nfinal (integer): Final tumour population size (default = 1000)
    - noise (string): Sequencing noise model (default = "betabinomial")
    - lower_cutoff (float): Lower frequency (range 0 - 0.5) bound for detecting subclones (default = 0.09)
    - upper_cutoff (float): Upper frequency (range 0 - 0.5) bound for detecting subclones (default = 0.41)
    - nsubclone_min (integer): Minimum number of subclones within (lower_cutoff, upper_cutoff) to accept for a positive selection simulation
    - nsubclone_max (integer): Maximum number of subclones within (lower_cutoff, upper_cutoff) to accept for a positive selection simulation
    - alt_reads (integer): Minimum number of alternate reads to call a mutation during virtual sequencing

Returns:
    - nested array: Positive selection simulation output
        - string: Simulation id
        - array: VAF, reads, depth for each mutation where each entry is an array with a length equal to number of detectable mutations
        - array: Subclone/driver fitness, timing, frequency information
        - array: Mutation count information
        - array: Simulation parameters
    - nested array: Neutral evolution simulation output
        - string: Simulation id
        - array: VAF, reads, depth for each mutation where each entry is an array with a length equal to number of detectable mutations
        - array: Subclone/driver fitness, timing, frequency information
        - array: Mutation count information
        - array: Simulation parameters
"""
function autoSimulation(b, d; Nfinal::Int64 = 1000, noise::String = "betabinomial", lower_cutoff::Float64 = 0.09, upper_cutoff::Float64 = 0.41, nsubclone_min::Int64 = 1, nsubclone_max::Int64 = 2, alt_reads::Int64 = 2)

    # Run positive selection simulation
    u, phi_b, nclonal, ndrivers, lambda, depth, rho = sampleParameters()
    println("{A} Starting positive simulation")
    println("--------------------------------")
    println("")
    n_detectable, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters = 0, [], [], [], []
    counter = 0
    while true
        counter += 1
        if counter == 10
            u, phi_b, nclonal, ndrivers, lambda, depth = sampleParameters()
            counter = 0
            println("Parameters failed to generate subclone number of interest, resampling parameters...")
        end
        n_detectable, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters = simulation(b, d, u, Nfinal = Nfinal, phi_b = phi_b, phi_d = 0.0, nclonal = nclonal, depth = depth, noise = noise, rho = rho, lower_cutoff = lower_cutoff, upper_cutoff = upper_cutoff, ndrivers = ndrivers, lambda = lambda, alt_reads = alt_reads)
        if nsubclone_max >= n_detectable >= nsubclone_min
            println("")
            println("Correct number of subclones generated for synthetic tumour with positive selection...")
            println("")
            break
        end
        println("Non-specified number of subclones generated, re-simulating positive selection...")
    end

    # Generate a synthetic neutral simulation using matched mutational burden parameters from positive selection simulation
    println("{B} Generating parameter paired fully synthetic neutral data")
    println("------------------------------")
    println("")

    scale_pos, scale_clonal = rand(Uniform(1, 1.5)), rand(Uniform(0.5, 1.5)) # Add variability to clonal and subclonal mutation counts in neutral synthetic tumours
    trim_tail = ifelse(rand(Uniform(0, 1), 1)[1] < 0.05, 1, 0) # Randomly remove tail from neutral distributions to inject more sparsity in neutral samples
    n_detectable, neutral_VAF, neutral_haplotypes, neutral_mutations, neutral_parameters = full_synthetic_neutral(depth = depth, Nfinal = Nfinal, nclonal = Int(round(nclonal * scale_clonal)), match_pos_muts = Int(round(positive_mutations[5] * scale_pos)), noise = noise, rho = rho, alt_reads = alt_reads, trim_tail = trim_tail)

    println("Saving synthetic data...")
    simid = randstring(10)
    positive_id, neutral_id = join(["P_", simid]), join(["N_", simid])

    return [positive_id, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters], [neutral_id, neutral_VAF, neutral_haplotypes, neutral_mutations, neutral_parameters]
end

"""
Generates a histogram of the mutation frequency information for a given tumour evolution simulation.

Input:
    - CanEvolve simulation output: A single simulation object e.g. p from p, n = autoSimulation(1, 0.1)
    - bins (integer): Number of bins in histogram

Returns:
    - Plots.Plot{Plots.GRBackend}: A histogram of the VAF information from a simulation. If sim_object is positive selection, blue lines indicate subclone frequency (cellular fraction / 2)
"""
function plotVAF(sim_object; bins::Int64 = 100)
    VAF = sim_object[2][1]
    scaler = fit(Histogram, VAF, nbins = bins)
    scaler = normalize(scaler, mode = :pdf)
    plot(scaler, fillcolor = RGB(100/255,180/255,137/255), linecolor = RGB(192/255,193/255,194/255), linealpha = 0.8, legend = false, xlims = (0,1), xtickfontsize=8,ytickfontsize=8)
    if sum(sim_object[3][1]) != 0
        vline!(sim_object[3][1], color = RGB(118/255,140/255,195/255), alpha = 1, linewidth = 3)
    end
    xlabel!("Variant allele frequency (VAF)")
    ylabel!("Density")
    density!(VAF, linecolor = RGBA(249/255,101/255,103/255,0), linewidth = 5, linealpha = 0.95)
end

"""
Computes frequency threshold under a binomial sequencing noise, where variance = np(1-p), model where the lower cutoff is variable to account for additional variability in empirical data

Input:
    - depth (integer/float): Mean sequencing depth
    - alt_reads (integer): Minimum number of reads covering mutation (default = 2)

Returns:
    - float: The lower frequency cutoff for VAFs that determines the limit of detection prior to adding sequencing noise.
    - float: The upper frequency cutoff that is fixed at 3 binomial standard deviations from 0.5. ## Note, this is generally not used in this framework
"""
function frequencyThresholds(depth; alt_reads::Int64 = 2)
    tail_shift = rand(Uniform(1, 3), 1)[1]
    f_min = (alt_reads/depth) + ( ( tail_shift*sqrt(alt_reads*(1-(alt_reads/depth))) ) / depth)
    f_max = 0.5 - ( ( 3.0*sqrt((0.5*depth)*(1-0.5)) ) / depth)
    return f_min, f_max
end

"""
Generates a histogram representation of the VAF distribution.

Input:
    - array: VAF information for each mutation
    - depth (integer): Mean sequencing depth (default = 100)
    - alt_reads (integer): Minimum number of reads covering mutation (default = 2)
    - k (integer): Number of bins to segment the VAF distribution
    - range (array): An array indicating lower and upper bound to constrain histogram range (default = [0.02, 0.5])
    - cutoff (integer): A 0 or 1 indicating if the lower bound should be based on range[1] or on 2 binomial standard devations above alt_reads/depth

Retuns:
    - array: An array of length k where each entry indicates the number of mutations in bins of width (range[2] - range[1])/k
"""
function uniformdensity(vaf; depth::Int64 = 100, alt_reads::Int64 = 2, k::Int64 = 100, range = [0.02, 0.5], cutoff::Int64 = 1)
    if cutoff == 1
        f_min = (alt_reads/depth) + ( ( 2*sqrt(alt_reads*(1-(alt_reads/depth))) ) / depth)
        h = np.histogram(vaf[(range[2] .> vaf .> f_min)], range = range, bins = k)
        nd = h[1]
    else
        h = np.histogram(vaf[(range[2] .> vaf .> range[1])], range = range, bins = k)
        nd = h[1]
    end
    return nd
end

"""
Normalizes a vector by dividing all entries by the maximum value in the vector.
"""
function rescale!(arr)
    mx = maximum(arr)
    arr = arr ./ mx
    return arr
end

"""
Transformation of the VAF information into normalized feature vectors.

Input:
    - CanEvolve simulation output: A single simulation object e.g. p from p, n = autoSimulation(1, 0.1)
    - lower_cutoff (float): Lower frequency (range 0 - 0.5) bound for detecting subclones (default = 0.09)
    - upper_cutoff (float): Upper frequency (range 0 - 0.5) bound for detecting subclones (default = 0.41)
    - k (array): Number of bins in histogram feature vector. Multiple entries of different sizes will generate multiple feature vectors (default = [100])

Returns:
    - array: Features
    - array: Labels (mode, nsubclones, subclone frequency(ies), subclone emergence time, subclone fitness, subclone fitness relative to population average, mutation rate, mean depth, sequencing overdispersion)
"""
function engineer(sim_object; lower_cutoff::Float64 = 0.09, upper_cutoff::Float64 = 0.41, k::Array{Int64} = [100])

    VAF, depth, alt_reads = sim_object[2][1], Int64(round(mean(sim_object[2][3]))), sim_object[5][end]

    features = []
    for bin_number in k
        ud1 = rescale!(uniformdensity(VAF, k = bin_number, alt_reads = alt_reads, range = [0.02, 0.5]))
        push!(features, ud1)
    end

    mode = ifelse(sum(sim_object[3][1]) == 0.0, 0, 1)
    if mode == 1
        detectable = intersect(findall(2*upper_cutoff .> sim_object[3][2] .> 2*lower_cutoff), findall(upper_cutoff .> sim_object[3][1] .> lower_cutoff))
        labels = [1, length(sim_object[3][1][detectable]), sim_object[3][1][detectable], sim_object[3][2][detectable], sim_object[3][4][detectable], sim_object[3][5][detectable], sim_object[5][3], depth, sim_object[5][10]]
    elseif mode == 0
        labels = [0, 0, 0, 0, 0, 0, 0, depth, sim_object[5][10]]
    else
        println("Did not identify proper mode.")
    end

    return features, labels
end

"""
Transformation of the VAF information into nonnormalized feature vectors.

Input:
    - CanEvolve simulation output: A single simulation object e.g. p from p, n = autoSimulation(1, 0.1)
    - lower_cutoff (float): Lower frequency (range 0 - 0.5) bound for detecting subclones (default = 0.09)
    - upper_cutoff (float): Upper frequency (range 0 - 0.5) bound for detecting subclones (default = 0.41)
    - k (array): Number of bins in histogram feature vector. Multiple entries of different sizes will generate multiple feature vectors (default = [100])

Returns:
    - array: Features
    - array: Labels (mode, nsubclones, subclone frequency(ies), subclone emergence time, subclone fitness, subclone fitness relative to population average, mutation rate, mean depth, sequencing overdispersion)
"""
function engineer_un(sim_object; lower_cutoff::Float64 = 0.09, upper_cutoff::Float64 = 0.41, k::Array{Int64} = [100])

    VAF, depth, alt_reads = sim_object[2][1], Int64(round(mean(sim_object[2][3]))), sim_object[5][end]
    dp = sim_object[5][8]
    features = []
    for bin_number in k
        ud1 = uniformdensity(VAF, k = bin_number, depth = dp, alt_reads = alt_reads, range = [0.02, 0.5], cutoff = 1)
        if length(k) == 1
            append!(features, ud1)
        else
            push!(features, ud1)
        end
    end

    mode = ifelse(sum(sim_object[3][1]) == 0.0, 0, 1)
    if mode == 1
        detectable = intersect(findall(2*upper_cutoff .> sim_object[3][2] .> 2*lower_cutoff), findall(upper_cutoff .> sim_object[3][1] .> lower_cutoff))
        labels = [1, length(sim_object[3][1][detectable]), sim_object[3][1][detectable], sim_object[3][3][detectable], sim_object[3][4][detectable], sim_object[3][5][detectable], sim_object[5][3], depth, sim_object[5][10]]
    elseif mode == 0
        labels = [0, 0, 0, 0, 0, 0, sim_object[5][3], depth, sim_object[5][10]]
    else
        println("Did not identify proper mode.")
    end

    return features, labels
end

"""
A replicate of autoSimulation but with fixed mutation rate sequencing depth and overdispersion (rho) parameters.

See autoSimulation for a full description of input parameters.
"""
function autoSimulation_fixed(b, d; u::Int64 = 20, depth::Int64 = 100, Nfinal::Int64 = 1000, noise::String = "betabinomial", rho::Float64 = 0.0, lower_cutoff::Float64 = 0.1, upper_cutoff::Float64 = 0.4, nsubclone_min::Int64 = 1, nsubclone_max::Int64 = 2, alt_reads::Int64 = 2)

    # Run simulation
    _, phi_b, nclonal, ndrivers, lambda, _ = sampleParameters()
    println("{A} Starting positive simulation")
    println("--------------------------------")
    println("")
    n_detectable, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters = 0, [], [], [], []
    counter = 0
    while true
        counter += 1
        if counter == 10
            _, phi_b, nclonal, ndrivers, lambda, _ = sampleParameters()
            counter = 0
            println("Parameters failed to generate subclone number of interest, resampling parameters...")
        end
        n_detectable, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters = simulation(b, d, u, Nfinal = Nfinal, phi_b = phi_b, phi_d = 0.0, nclonal = nclonal, depth = depth, noise = noise, rho = rho, lower_cutoff = lower_cutoff, upper_cutoff = upper_cutoff, ndrivers = ndrivers, lambda = lambda, alt_reads = alt_reads)
        if nsubclone_max >= n_detectable >= nsubclone_min
            println("")
            println("Correct number of subclones generated for synthetic tumour with positive selection...")
            println("")
            break
        end
        println("Non-specified number of subclones generated, re-simulating positive selection...")
    end

    println("{B} Generating parameter paired fully synthetic neutral data")
    println("------------------------------")
    println("")

    scale_pos, scale_clonal = rand(Uniform(1, 1.5)), rand(Uniform(0.5, 1.5)) # Add variability to clonal and subclonal mutation counts in neutral synthetic tumours
    trim_tail = ifelse(rand(Uniform(0, 1), 1)[1] < 0.05, 1, 0) # Randomly remove tail from neutral distributions to inject more sparsity in neutral samples
    n_detectable, neutral_VAF, neutral_haplotypes, neutral_mutations, neutral_parameters = full_synthetic_neutral(depth = depth, Nfinal = Nfinal, nclonal = Int(round(nclonal * scale_clonal)), match_pos_muts = Int(round(positive_mutations[5] * scale_pos)), noise = noise, rho = rho, alt_reads = alt_reads, trim_tail = trim_tail)

    println("Saving synthetic data...")
    # Add identifiers for future mapping if needed
    simid = randstring(10)
    positive_id, neutral_id = join(["P_", simid]), join(["N_", simid])

    return [positive_id, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters], [neutral_id, neutral_VAF, neutral_haplotypes, neutral_mutations, neutral_parameters]
end

"""
A replicate of autoSimulation but with fixed mutation rate and sequencing depth parameters.

See autoSimulation for a full description of input parameters.
"""
function autoSimulation_fixed_depth(b, d; depth::Int64 = 100, Nfinal::Int64 = 1000, noise::String = "betabinomial", lower_cutoff::Float64 = 0.1, upper_cutoff::Float64 = 0.4, nsubclone_min::Int64 = 1, nsubclone_max::Int64 = 2, alt_reads::Int64 = 2)

    # Run simulation
    u, phi_b, nclonal, ndrivers, lambda, _, rho = sampleParameters()
    println("{A} Starting positive simulation")
    println("--------------------------------")
    println("")
    n_detectable, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters = 0, [], [], [], []
    counter = 0
    while true
        counter += 1
        if counter == 10
            _, phi_b, nclonal, ndrivers, lambda, _ = sampleParameters()
            counter = 0
            println("Parameters failed to generate subclone number of interest, resampling parameters...")
        end
        n_detectable, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters = simulation(b, d, u, Nfinal = Nfinal, phi_b = phi_b, phi_d = 0.0, nclonal = nclonal, depth = depth, noise = noise, rho = rho, lower_cutoff = lower_cutoff, upper_cutoff = upper_cutoff, ndrivers = ndrivers, lambda = lambda, alt_reads = alt_reads)
        if nsubclone_max >= n_detectable >= nsubclone_min
            println("")
            println("Correct number of subclones generated for synthetic tumour with positive selection...")
            println("")
            break
        end
        println("Non-specified number of subclones generated, re-simulating positive selection...")
    end

    println("{B} Generating parameter paired fully synthetic neutral data")
    println("------------------------------")
    println("")

    scale_pos, scale_clonal = rand(Uniform(1, 1.5)), rand(Uniform(0.5, 1.5)) # Add variability to clonal and subclonal mutation counts in neutral synthetic tumours
    trim_tail = ifelse(rand(Uniform(0, 1), 1)[1] < 0.05, 1, 0) # Randomly remove tail from neutral distributions to inject more sparsity in neutral samples
    n_detectable, neutral_VAF, neutral_haplotypes, neutral_mutations, neutral_parameters = full_synthetic_neutral(depth = depth, Nfinal = Nfinal, nclonal = Int(round(nclonal * scale_clonal)), match_pos_muts = Int(round(positive_mutations[5] * scale_pos)), noise = noise, rho = rho, alt_reads = alt_reads, trim_tail = trim_tail)

    println("Saving synthetic data...")
    # Add identifiers for future mapping if needed
    simid = randstring(10)
    positive_id, neutral_id = join(["P_", simid]), join(["N_", simid])

    return [positive_id, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters], [neutral_id, neutral_VAF, neutral_haplotypes, neutral_mutations, neutral_parameters]
end
