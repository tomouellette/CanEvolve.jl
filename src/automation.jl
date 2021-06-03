# ---------------------------------------------------------------------------
# automation.jl: Various functions for synthetic data generation and analysis
# ---------------------------------------------------------------------------
np = pyimport("numpy")

# Note: manual specification of parameter ranges
function sampleParameters()
    u = rand(DiscreteUniform(10, 500), 1)[1]
    phi_b = rand(Uniform(0.001, 0.1), 1)[1] / u # Normalize by u so time of driver occurences are similar across different mutation rates
    nclonal = rand(DiscreteUniform(100, u*10), 1)[1] # We remove clonal during inference so just keep scaled relative to mutation rate
    ndrivers = rand(DiscreteUniform(1, 3), 1)[1]
    lambda = rand(Uniform(1, 3), 1)[1]
    depth = rand(DiscreteUniform(50, 125), 1)[1]
    return u, phi_b, nclonal, ndrivers, lambda, depth
end

function autoSimulation(b, d; Nfinal::Int64 = 10000, noise::String = "binomial", rho::Float64 = 0.0, lower_cutoff::Float64 = 0.1, upper_cutoff::Float64 = 0.4, nsubclone_min::Int64 = 1, nsubclone_max::Int64 = 2, alt_reads::Int64 = 2)

    # Run simulation
    u, phi_b, nclonal, ndrivers, lambda, depth = sampleParameters()
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

    println("{B} Generating parameter paired fully synthetic neutral data")
    println("------------------------------")
    println("")

    n_detectable, neutral_VAF, neutral_haplotypes, neutral_mutations, neutral_parameters = full_synthetic_neutral(depth = depth, Nfinal = Nfinal, nclonal = nclonal, match_pos_muts = positive_mutations[5], noise = noise, rho = rho, alt_reads = alt_reads)

    println("Saving synthetic data...")
    # Add identifiers for future mapping if needed
    simid = randstring(10)
    positive_id, neutral_id = join(["P_", simid]), join(["N_", simid])

    return [positive_id, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters], [neutral_id, neutral_VAF, neutral_haplotypes, neutral_mutations, neutral_parameters]
end

function plotVAF(sim_object; bins::Int64 = 100)
    VAF = sim_object[2][1]
    scaler = fit(Histogram, VAF, nbins = bins)
    scaler = normalize(scaler, mode = :pdf)
    plot(scaler, fillcolor = "black", linecolor = "gray", linealpha = 1, legend = false, xlims = (0,1), xtickfontsize=8,ytickfontsize=8)
    if sum(sim_object[3][1]) != 0
        vline!(sim_object[3][1], color = "blue", alpha = 0.8, linewidth = 3)
    end
    xlabel!("Variant allele frequency (VAF)")
    ylabel!("Density")
    density!(VAF, linecolor = "red", linewidth = 5, linealpha = 0.8)
end

function frequencyThresholds(depth; alt_reads::Int64 = 2)
    f_min = (alt_reads/depth) + ( ( 2.0*sqrt(alt_reads*(1-(alt_reads/depth))) ) / depth)
    f_max = 0.5 - ( ( 3.0*sqrt((0.5*depth)*(1-0.5)) ) / depth)
    return f_min, f_max
end

function mapfrequency(udensity; depth::Int64 = 100, alt_reads::Int64 = 2, include_clonal = false)
    f_min, f_max = frequencyThresholds(depth, alt_reads = alt_reads)
    if include_clonal == true
        width = (0.5 - f_min) / length(udensity)
    else
        width = (f_max - f_min) / length(udensity)
    end
    freqmap = [f_min + width*i for i in 1:length(udensity)]
    return freqmap
end

function uniformdensity(vaf; k::Int64 = 100, depth::Int64 = 50, alt_reads::Int64 = 2, include_clonal = false)
    f_min, f_max = frequencyThresholds(depth, alt_reads = alt_reads)
    if include_clonal == true
        h = np.histogram(vaf[(0.5 .> vaf .> f_min)], bins = k, density = true)
    else
        h = np.histogram(vaf[(f_max .> vaf .> f_min)], bins = k, density = true)
    end
    nd = h[1] ./ np.sum(h[1])
    return nd
end

function rescale(arr)
    mx = maximum(arr)
    arr = arr ./ mx
    return arr
end

function engineer(sim_object; lower_cutoff::Float64 = 0.09, upper_cutoff::Float64 = 0.41, k::Array{Int64} = [100])

    VAF, depth, alt_reads = sim_object[2][1], Int64(round(mean(sim_object[2][3]))), sim_object[5][end]

    features = []
    for bin_number in k
        ud1 = rescale(uniformdensity(VAF, k = bin_number, depth = depth, alt_reads = alt_reads, include_clonal = false))
        mf1 = 2 .* mapfrequency(ud1, depth = depth, alt_reads = alt_reads, include_clonal = false)
        ud2 = rescale(uniformdensity(VAF, k = bin_number, depth = depth, alt_reads = alt_reads, include_clonal = true))
        mf2 = 2 .* mapfrequency(ud2, depth = depth, alt_reads = alt_reads, include_clonal = true)
        feature = transpose(hcat(ud1, mf1, ud2, mf2))
        push!(features, feature)
    end

    # Labels: MODE (0, 1), NSUBCLONES (0, 1, 2), SC_FREQ, SC_FITNESS, SC_RELATIVE_AVG, MUTRATE, DEPTH, RHO
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

# Additional functions
function autoSimulation_fixed(b, d; u::Int64 = 20, depth::Int64 = 100, Nfinal::Int64 = 10000, noise::String = "binomial", rho::Float64 = 0.0, lower_cutoff::Float64 = 0.1, upper_cutoff::Float64 = 0.4, nsubclone_min::Int64 = 1, nsubclone_max::Int64 = 2, alt_reads::Int64 = 2)

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

    n_detectable, neutral_VAF, neutral_haplotypes, neutral_mutations, neutral_parameters = full_synthetic_neutral(depth = depth, Nfinal = Nfinal, nclonal = nclonal, match_pos_muts = positive_mutations[5], noise = noise, rho = rho, alt_reads = alt_reads)

    println("Saving synthetic data...")
    # Add identifiers for future mapping if needed
    simid = randstring(10)
    positive_id, neutral_id = join(["P_", simid]), join(["N_", simid])

    return [positive_id, positive_VAF, positive_haplotypes, positive_mutations, positive_parameters], [neutral_id, neutral_VAF, neutral_haplotypes, neutral_mutations, neutral_parameters]
end
