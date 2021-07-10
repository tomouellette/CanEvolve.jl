# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# simulation.jl: Tumour evolution under a stochastic branching process
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
This simulation framework is adapted/extended from CancerSeqSim.jl by Marc Williams (Williams et al. 2018), as well as ideas from the rkMC in Waclaw et al. 2015.

Our framework differs from previous frameworks by:
   (i) Allowing for stochastic arrival of driver, and/or passenger, mutations that scale growth rates (b - d)
   (i) Integrating a multiplicative or additive fitness model that enables an arbitrary number of unique driver haplotypes to arise (similar to McFarland et al. 2014)
   (iii) Allow for deleterious passenger mutations that decrease fitness if interested
   (iv) At small N, we use a fully synthetic generation of neutral evolution data as genetic drift can lead to excess spurious subclones that lead to a biased "ground truth" dataset for ML inferences
"""

# =========================================================================================================================================================================================================

"""
A class to store each cell's mutation and evolutionary information.
"""
mutable struct cell

    neutral::Array{Int64,1}
    beneficial::Array{Int64,1}
    deleterious::Array{Int64,1}
    s_b::Array{Float64,1}
    s_d::Array{Float64,1}
    time_of_birth::Float64

end

"""
Takes in a cell and adds mutations.

Input:
    - struct: An object of struct cell
    - integer: Mutation rate (per genome per division)
    - float: Probability of subclone/driver event (per division probability)
    - float: Probability of deleterious mutation (per division probability)
    - integer: A mutation id to enable counting post-simulation
    - ndrivers: The current number of subclone/driver events remaining in the simulation
    - lambda (float): Mean of exponential selection coefficient distribution (default = 0.1)
    - fitness_model (string): Increase fitness under an "additive" or "multiplicative" model (default = "multiplicative")

Returns:
    - struct: An updated object of struct cell
    - integer: An updated mutation id that increased by the number of new mutations generated
    - float: The fitness scaled birth rate of the cell, > 1 if a new driver event occured else 1
    - float: The fitness scaled death rate of the cell, < 1 if a new deleterious passenger occured else 1
    - float: The selection coefficients of the new mutations
    - integer: Number of driver subclone/driver events remaining in the simulation

"""
function mutate(cell, µ, phi_b, phi_d, mut_id, ndrivers, lambda = 0.1, fitness_model = "multiplicative")

    neutral = rand(Poisson(µ * (1 - phi_b - phi_d)))
    beneficial = rand(Poisson(µ * phi_b))
    deleterious = rand(Poisson(µ * phi_d))

    fitness = [0.0 for i in 1:neutral]

    cell.neutral = append!(cell.neutral, mut_id:mut_id+neutral-1)
    mut_id = mut_id + neutral

    update_bmax, update_dmax = 1.0, 1.0
    if (beneficial > 0) & (ndrivers > 0)
        cell.beneficial = append!(cell.beneficial, mut_id:mut_id+beneficial-1)
        sb = rand(Exponential(lambda), beneficial)
        cell.s_b, fitness = append!(cell.s_b, sb), append!(fitness, sb)
        mut_id = mut_id + beneficial
        if fitness_model == "multiplicative"
            update_bmax = prod(1 .+ cell.s_b)
        elseif fitness_model == "additive"
            update_bmax = sum(1 .+ cell.s_b)
        else
            update_bmax = prod(1 .+ cell.s_b)
        end
        ndrivers = ndrivers - beneficial
    end

    if deleterious > 0
        cell.deleterious = append!(cell.deleterious, mut_id:mut_id+deleterious-1)
        sd = -0.001 .- rand(Exponential(0.001), deleterious)
        cell.s_d, fitness = append!(cell.s_d, sd), append!(fitness, sd)
        mut_id = mut_id + deleterious
        update_dmax = prod(1 .+ cell.s_d)
        if fitness_model == "multiplicative"
            update_dmax = prod(1 .+ cell.s_d)
        elseif fitness_model == "additive"
            update_dmax = sum(1 .+ cell.s_d)
        else
            update_dmax = prod(1 .+ cell.s_d)
        end
    end

    return cell, mut_id, update_bmax, update_dmax, fitness, ndrivers

end

"""
Primes the simulation with values and arrays for tracking events.
"""
function primeSim(;nclonal::Int64 = 1)

    # One cell at time 0
    N, t, mut_id, = 1, 0.0, nclonal + 1

    # Keep track of time, population size, and cells
    track_N = Int64[]
    track_t = Float64[]
    track_mut_time = Float64[]
    track_cells = cell[]
    track_fitness = Float64[]

    # Add initial values to arrays
    push!(track_N, N)
    push!(track_t, t)
    append!(track_mut_time, [t for i in 1:nclonal])
    append!(track_fitness, [0.0 for i in 1:nclonal])
    if nclonal > 0
        push!(track_cells, cell([i for i in 1:nclonal],[],[],[],[],t))
    else
        push!(track_cells, cell([],[],[],[],[],t))
    end

    return N, track_N, t, track_t, track_mut_time, track_cells, mut_id, track_fitness

end

"""
See autoSimulation documentation in automation.jl for description of parameters.
"""
function evolveTumour(b, d, µ; Nfinal::Int64 = 10000, phi_b::Float64 = 0.01, phi_d::Float64 = 0.5, nclonal::Int64 = 0, ndrivers::Int64 = 50, lambda::Float64 = 0.1, fitness_model::String = "multiplicative")

    bd_max = b + d
    update_bmax, update_dmax = 1, 1
    N, track_N, t, track_t, track_mut_time, track_cells, mut_id, track_fitness = primeSim(nclonal = nclonal)

    ncounter = 1
    while N < Nfinal

        selectcell = rand(1:N)
        r = rand(Uniform(0, bd_max))

        cell_fitness = prod(1 .+ track_cells[selectcell].s_b) * prod(1 .+ track_cells[selectcell].s_d)
        birthrate_cell = cell_fitness*b
        deathrate_cell = cell_fitness*d
        mut_id_t = mut_id # Store this id to add time of occurence each mutation

        if birthrate_cell > r

            # Cell undergoes division C -> 2C
            N = N + 1
            push!(track_N, N)
            push!(track_t, t)
            push!(track_cells, deepcopy(track_cells[selectcell]))

            # Add new mutations to both cells
            track_cells[selectcell], mut_id, ubm1, udm1, add_fitness1, ndrivers = mutate(track_cells[selectcell], µ, phi_b, phi_d, mut_id, ndrivers, lambda, "multiplicative")
            track_cells[end], mut_id, ubm2, udm2, add_fitness2, ndrivers = mutate(track_cells[end], µ, phi_b, phi_d, mut_id, ndrivers, lambda, "multiplicative")

            append!(track_fitness, vcat(add_fitness1, add_fitness2))
            append!(track_mut_time, repeat(t:t, mut_id - mut_id_t))
            track_cells[selectcell].time_of_birth = t
            track_cells[end].time_of_birth = t

            #Δt =  - 1/(bd_max * (N-1)) * log(rand())
            Δt = 1/(N-1)
            t = t + Δt

            #update_bmax, update_dmax = max(ubm1, ubm2, update_bmax), min(ubm1, ubm2, update_dmax)
            update_bmax, update_dmax = max(ubm1, ubm2, update_bmax), min(udm1, udm2, update_dmax)
            bd_max = update_bmax*update_dmax*b + update_bmax*update_dmax*d

        end

        if birthrate_cell + deathrate_cell > r >= birthrate_cell

            # Cell dies
            N = N - 1

            if N == 0
                println("Restarting simulation. Death rate must be too high.")
                bd_max = b + d
                update_bmax, update_dmax = 1, 1
                N, track_N, t, track_t, track_mut_time, track_cells, mut_id, track_fitness = primeSim(nclonal = nclonal)
            else
                deleteat!(track_cells, selectcell)
                push!(track_N, N)
                push!(track_t, t)

                #Δt = - 1/(bd_max * (N+1)) * log(rand())
                Δt = 1/(N+1)
                t = t + Δt
            end

        end

        # If r equals the maximum, nothing happens
        if birthrate_cell + deathrate_cell <= r

            push!(track_N, N)
            #Δt =  -1/(bd_max * N) * log(rand())
            Δt = 1/N
            t = t + Δt
            push!(track_t, t)

        end

    end

    return N, track_N, t, track_t, track_mut_time ./ t, track_cells, (mut_id - 1), track_fitness # Return ID - 1 as mut_id is primed for the next iteration (therefore it's n + 1)

end

"""
Beta binomial distribution (reference: CancerSeqSim.jl by Marc Williams)
"""
function BetaBin(n, p, rho)
    # Note: P(X = k) = {n \\choose k} B(k + \\alpha, n - k + \\beta) / B(\\alpha, \\beta),  \\quad \\text{ for } k = 0,1,2, \\ldots, n.
    mu = p * n
    alpha = (mu / n) * ((1 / rho) - 1)
    beta = n * alpha / mu - alpha
    return rand(Binomial(n, rand(Beta(alpha, beta))))
end

"""
Virtually sequences a tumour where depth is binomially distributed and reads are binomial or beta-binomially distributed.

Input:
    - integer: Mean sequencing depth specified for simulation
    - integer: Final population size of simulated tumour
    - array: VAF for each mutation
    - array: Fitness/selection coefficient of each mutation
    - array: Functional type for each mutation e.g. "N", "NS", "S" for neutral, nonsynonymous, or synonymous (##Note: not generally used in this framework)
    - limit_of_detection (float): The lower VAF bound detectable prior to adding sequencing noise.
    - noise (string): Sequencing noise model either "binomial" or "betabinomial" (default = "betabinomial")
    - rho (float): Overdispersion parameter for betabinomial distribution
    - alt_reads (integer): Minimum number of alternate reads to call a mutation during virtual sequencing

Returns:
    - array: VAF for each mutation injected with sequencing noise
    - array: Reads for each mutation
    - array: Sequencing depth for each mutation
    - array: Fitness of each mutation where mutations that were lost during virtual sequencing are removed
    - array: Functional type for each mutation where mutations that were lost during virtual sequencing are removed
"""
function sequence(depth, N, VAF, fitness, mut_type; limit_of_detection::Float64 = 0.1, noise::String = "betabinomial", rho::Float64 = 0.0025, alt_reads::Int64 = 2)

    # Binomially distributed depth with binomially distributed reads
    VAF = VAF ./ 2 # Takes cell fractions and converts to frequency
    VAF, fitness, mut_type = VAF[VAF .> limit_of_detection], fitness[VAF .> limit_of_detection], mut_type[VAF .> limit_of_detection]
    dp = rand.(Binomial.(N, repeat(depth:depth, length(VAF)) ./ N))
    if (noise == "binomial") | (rho == 0.0)
        reads = rand.(Binomial.(dp, VAF))
    end
    if (noise == "betabinomial") & (rho > 0.0)
        reads = BetaBin.(dp, VAF, rho)
    end

    vafadjust = reads ./ dp
    vafadjust = replace(vafadjust, NaN=>0)

    # Remove mutations below hard alternate read cutoff
    remove_zeros = findall(vafadjust .> (alt_reads / depth))

    return vafadjust[remove_zeros], reads[remove_zeros], dp[remove_zeros], fitness[remove_zeros], mut_type[remove_zeros]
end

"""
Parses information on detectable subclones/drivers.

Input:
    - array: Cells generated during virtual tumour growth and evolution
    - array: VAF for each mutation
    - array: Mutation ids that one-to-one map with VAF, mutation time, and mutation fitness information
    - integer: Final tumour population size
    - array: Time of emergence for each mutation
    - array: Fitness/selection coefficient of each mutation
    - lower_cutoff (float): Lower frequency (range 0 - 0.5) bound for detecting subclones (default = 0.09)
    - upper_cutoff (float): Upper frequency (range 0 - 0.5) bound for detecting subclones (default = 0.41)

Returns:
    - array: Subclone/haplotype ids for each detectable haplotype where an id is of the form [integer, integer, ...] where each integer represents the mutation id
    - array: Frequency of each haplotype (cellular fraction / 2)
    - array: Emergence times for each haplotype based on age of the most recent mutation in the haplotype
    - array: Multiplicative fitness of the haplotype
    - array: Fitness of the haplotype divided by the average population fitness
    - float: Average population fitness
    - float: Fitness of all cells that are not subclonal/driver haplotypes
    - array: Cellular fraction/proportions of each subclone/haplotype
    - integer: Number of detectable subclones i.e. subclone in (lower_cutoff, upper_cutoff)
    - integer: Number of undetectable subclones i.e. subclone in (0, lower_cutoff)U(upper_cutoff, 0.5)
"""
function detectableDrivers(track_cells, VAF, mut_ids, N, track_mut_time, fitness; lower_cutoff::Float64 = 0.1, upper_cutoff::Float64 = 0.4)

    if length(findall(fitness .> 0)) > 0
        selected = unique([i.beneficial for i in track_cells])
        selected = selected[isempty.(selected) .== false]
        ids = mut_ids[intersect(findall(VAF .> 2*lower_cutoff), findall(fitness .> 0))]
        detectable = []
        for i in 1:length(selected)
            v = selected[i]
            if all([j in ids for j in v]) == true
                push!(detectable, i)
            else
                continue
            end
        end

        # If no detectable driver haplotypes then simply return background fitness of entire population
        if length(detectable) < 1
            bgfitness = []
            for cell in track_cells
                push!(bgfitness, prod(1.0 .+ cell.s_b)*prod(1.0 .+ cell.s_d))
            end
            return [], [], [], [], [], [], [], mean(bgfitness), [], 0, 0
        end

        # Compute driver haplotype frequency
        detectable = selected[detectable]
        haplotype_counts = [0 for i in 1:length(detectable)]
        mfitness = [0.0 for i in 1:length(detectable)]
        bgfitness, bg = [], 0
        for cell in track_cells
            for i in 1:length(detectable)
                if all([j in cell.beneficial for j in detectable[i]]) == true
                    haplotype_counts[i] += 1
                    mfitness[i] += prod(1 .+ cell.s_b) * prod(1 .+ cell.s_d)
                else
                    bg = 1
                end
            end
            if bg == 1
                push!(bgfitness, prod(1 .+ cell.s_b) * prod(1 .+ cell.s_d))
                bg = 0
            end
        end

        mfitness = mfitness ./ haplotype_counts
        bgfitness = Any[]
        if length(bgfitness) > 0
            bgfitness = mean(bgfitness)
        else
            bgfitness = 1
        end

        # Find emergence times for each unique driver haplotypes
        times = [] # Find time of haplotype emergence based on most recent mutation in haplotype
        for haplotype in detectable
            ages = findmax(track_mut_time[collect(Iterators.flatten(haplotype))])[1]
            push!(times, ages)
        end

        N = 2*N

        # Compute cellular proportion from haplotype frequencies
        p = sortperm(length.(detectable)) # Sort all arrays by haplotype length
        detectable, haplotype_counts, haplotype_frequency, times, mfitness, bgfitness = detectable[p], haplotype_counts[p], haplotype_counts[p] ./ N, times[p], mfitness[p], bgfitness
        cellular_proportion = 2 .* adjustProportion(haplotype_frequency, detectable, mfitness, threshold = lower_cutoff)

        # Compute relative fitness based on average population fitness
        average_fitness = sum(cellular_proportion .* mfitness) + (bgfitness * (1 - sum(cellular_proportion)) )
        relative_fitness = mfitness ./ average_fitness

        # Filter simulations that have haplotypes in a detectable frequency range
        n_undetectable = length(detectable)
        detect_ind = intersect(findall(2*upper_cutoff .> cellular_proportion .> 2*lower_cutoff), findall(upper_cutoff .> haplotype_frequency .> lower_cutoff))
        #detectable, haplotype_counts, haplotype_frequency, times, mfitness, relative_fitness, cellular_proportion = detectable[detect_ind], haplotype_counts[detect_ind], haplotype_frequency[detect_ind], times[detect_ind], mfitness[detect_ind], relative_fitness[detect_ind], cellular_proportion[detect_ind]
        detectable_ = detectable[detect_ind]
        n_detectable = length(detectable_)

        return detectable, haplotype_counts, haplotype_frequency, times, mfitness, relative_fitness, average_fitness, bgfitness, cellular_proportion, n_detectable, n_undetectable

    else

        bgfitness = []
        for cell in track_cells
            push!(bgfitness, prod(1.0 .+ cell.s_b)*prod(1.0 .+ cell.s_d))
        end

        return [], [], [], [], [], [], [], mean(bgfitness), [], 0, 0

    end

end

"""
Get mutation information for each simulated tumour population.
"""
function processTumour(N, track_N, t, track_t, track_mut_time, track_cells, track_fitness, mut_id)

    # Get counts for each mutation in tumour population
    counts_neutral =  try sort(collect(countmap(reduce(vcat, [cell.neutral for cell in track_cells])))) catch e 0 end
    counts_beneficial = try sort(collect(countmap(reduce(vcat, [cell.beneficial for cell in track_cells])))) catch e 0 end
    counts_deleterious = try sort(collect(countmap(reduce(vcat, [cell.deleterious for cell in track_cells])))) catch e 0 end

    # Extract mutation ids and counts for each mutation type
    function aggregate(cntmap, name)
        if cntmap != 0
            ids, cnts, names = [i[1] for i in cntmap], [i[2] for i in cntmap], [name for i in 1:length(cntmap)]
        else
            ids, cnts, names = [], [], []
        end
        return ids, cnts, names
    end

    n_ids, n_counts, n_names = aggregate(counts_neutral, "N")
    b_ids, b_counts, b_names = aggregate(counts_beneficial, "B")
    d_ids, d_counts, d_names = aggregate(counts_deleterious, "D")

    # Compile data table for mutations, counts, and fitness
    mut_ids = vcat(n_ids, b_ids, d_ids)
    counts = vcat(n_counts, b_counts, d_counts)
    mut_type = vcat(n_names, b_names, d_names)
    VAF = vcat(n_counts ./ N, b_counts ./ N, d_counts ./ N)

    return mut_ids, track_mut_time[mut_ids], counts, mut_type, VAF, track_fitness[mut_ids]

end

"""
Adjust the subclone proportions (true cellular fractions) to account for nesting in the >1 subclone case
"""
function adjustProportion(dhap_freq, dhap_ids, mfitness; threshold::Float64 = 0.05)

    # Adjust haplotype frequencies to true cellular proportions (relevant if there are nested subclones)
    sort_ind = reverse(sortperm(length.(dhap_ids)))
    freqs = dhap_freq[sort_ind]
    haps = dhap_ids[sort_ind]
    mfitness = mfitness[sort_ind]

    update_freqs = [0.0 for i in 1:length(freqs)]
    for i in 1:length(haps)
        hap = haps[i]
        f = freqs[i]
        for j in 1:i
            if hap == haps[j]
                continue
            elseif all([k in haps[j] for k in hap])
                if update_freqs == 0.0
                    f = f - freqs[j]
                else
                    f = f - update_freqs[j]
                end
            else
                continue
            end
        end
        update_freqs[i] = f
    end

    return reverse(update_freqs)
end

function mutationMetrics(VAF, mut_type, fitness; f_min::Float64 = 0.1, f_max::Float64 = 0.4)

    beneficial = sum(fitness .> 0.0)
    deleterious = sum(fitness .< 0.0)
    neutral = sum(mut_type .== "N")
    total = length(VAF)
    subclonal = length(VAF[VAF .< f_max])

    return beneficial, deleterious, neutral, total, subclonal

end

"""
See autoSimulation documentation in automation.jl for description of parameters.
"""
function simulation(b, d, u; Nfinal::Int64 = 10000, phi_b::Float64 = 0.0, phi_d::Float64 = 0.0, nclonal::Int64 = 100, depth::Int64 = 100, noise::String = "betabinomial", rho::Float64 = 0.0, lower_cutoff::Float64 = 0.1, upper_cutoff::Float64 = 0.4, ndrivers::Int64 = 20, lambda::Float64 = 0.1, fitness_model::String = "multiplicative", alt_reads::Int64 = 2)

    f_min, f_max = frequencyThresholds(depth, alt_reads = alt_reads)

    println("[1] Growing tumour...")
    @time N, track_N, t, track_t, track_mut_time, track_cells, mut_id, track_fitness = evolveTumour(b, d, u, Nfinal = Nfinal, phi_b = phi_b, phi_d = phi_d, nclonal = nclonal, ndrivers = ndrivers, lambda = lambda, fitness_model = fitness_model)
    println("")

    println("[2] Counting mutations...")
    @time mut_ids, mut_time, counts, mut_type, VAF, fitness = processTumour(N, track_N, t, track_t, track_mut_time, track_cells, track_fitness, mut_id)
    println("")

    println("[3] Finding detectable haplotypes/subclones...")
    if (phi_b == 0.0)
        @time haplotype_ids, haplotype_counts, haplotype_frequency, haplotype_times, absolute_fitness, relative_fitness, average_fitness, background_fitness, cellular_proportion, n_detectable, n_undetectable = [], [], [], [], [], [], 0.0, 0.0, [], 0, 0
    else
        @time haplotype_ids, haplotype_counts, haplotype_frequency, haplotype_times, absolute_fitness, relative_fitness, average_fitness, background_fitness, cellular_proportion, n_detectable, n_undetectable = detectableDrivers(track_cells, VAF, mut_ids, N, track_mut_time, fitness, lower_cutoff = lower_cutoff, upper_cutoff = upper_cutoff)
    end
    println("")

    println("[4] Virtually sequencing the tumour with an injection of ", noise, " sequencing noise...")
    @time VAF, reads, DP, fitness, mut_type = sequence(Int64(depth), Nfinal, VAF, fitness, mut_type, noise = noise, limit_of_detection = f_min, rho = rho, alt_reads = alt_reads)
    println("")

    println("[5] Getting mutation info for saving...")
    @time beneficial, deleterious, neutral, total, subclonal = mutationMetrics(VAF, mut_type, fitness, f_min = f_min, f_max = f_max)
    println("")

    println("[6] Saving parameters...")
    parameters = [b, d, u, Nfinal, phi_b, phi_d, nclonal, depth, noise, rho, lower_cutoff, upper_cutoff, ndrivers, lambda, fitness_model, alt_reads]

    return n_detectable, [VAF, reads, DP], [haplotype_frequency, cellular_proportion, haplotype_times, absolute_fitness, relative_fitness], [beneficial, deleterious, neutral, total, subclonal], parameters
end

"""
Generates neutral evolution synthetic VAF data by sampling from Pareto distributions with empirically realistic shape and and scale parameters.

Notes:
    - The Pareto distribution that defines the neutral VAF tail is parameterized by shape and scale
    - Caravagna et al. 2020 fit neutral tails, i.e. a Pareto distribution, to many PCAWG samples.
    - To generate realistic shape and scale values during synthetic neutral data generaton,
      we built sampling distributions for the shape and scale parameters by performing MLE fits
      to all PCAWG samples above 50x coverage. (Although default scale parameter used here is f_min)

Input:
    - depth (integer): Mean sequencing depth (default = 100)
    - Nfinal (integer): Final tumour population size (default = 1000)
    - nclonal (integer): Number of clonal mutations  (default = 500)
    - match_pos_muts (integer): If pairing with a positive selection simulation, specifies number of non-clonal mutations to add (Default = 0; randomly sample mutation counts)
    - noise (string): Sequencing noise model either "binomial" or "betabinomial" (default = "betabinomial")
    - rho (float): Overdispersion parameter for betabinomial distribution
    - alt_reads (integer): Minimum number of alternate reads to call a mutation during virtual sequencing
    - trim_tail (integer): Specifies if tail of distribution should be randomly trimmed at 0.1 - 0.3 VAF. No trim if 0 or 1 to trim (default = 0)

Returns:
    - CanEvolve simulation object: A neutral evolution simulation object. See autoSimulation description in automation.jl  for more information
"""
function full_synthetic_neutral(;depth::Int64 = 100, Nfinal::Int64 = 1000, nclonal::Int64 = 500, match_pos_muts = 0, noise::String = "betabinomial", rho::Float64 = 0.002, alt_reads::Int64 = 2, trim_tail::Int64 = 0)

    """
    """
    # Fit of Caravagna et al 2020 shape parameters with gamma + E distribution: shape 4.016782, rate 6.404820
    shape = rand(Gamma(4.016782, 1/6.404820), 1)[1] + 0.6150006 # Scale = 1/rate

    # Fit of Caravagna et al. 2020 scale parameters with Pareto distribution: scale 0.05010 shape 28.43944
    #scale = rand(Pareto(28.43944, 0.05010), 1)[1] # We just use f_min for scale parameter in general

    # Match the number of subclonal mutations relative to matched positive stochastic simulation
    if match_pos_muts == 0
        nmuts = rand(DiscreteUniform(100, 20000), 1)[1]
    else
        nmuts = match_pos_muts
    end

    # Generate neutral tail VAFs
    f_min, f_max = frequencyThresholds(depth, alt_reads = alt_reads)
    vaf_nonclonal = rand(Pareto(shape, f_min), nmuts)
    vaf_nonclonal = vaf_nonclonal[findall(vaf_nonclonal .< 1)]
    vaf_clonal = [0.5 for i in 1:nclonal]
    vafs = np.hstack([vaf_nonclonal, vaf_clonal])

    # Synthetic sequencing
    if trim_tail == 1
        f_min = rand(Uniform(0.1, 0.3), 1)[1]
        vafs = vafs[vafs .> f_min]
    else
        vafs = vafs[vafs .> f_min]
    end
    dp = rand.(Binomial.(Nfinal, repeat(depth:depth, length(vafs)) ./ Nfinal))
    if (noise == "binomial") | (rho == 0.0)
        reads = rand.(Binomial.(dp, vafs))
    end
    if (noise == "betabinomial") & (rho > 0.0)
        reads = BetaBin.(dp, vafs, rho)
    end
    vafs = reads ./ dp
    vafs = replace(vafs, NaN=>0)

    # Remove mutations below hard alternate read cutoff
    remove_zeros = findall(vafs .> (alt_reads / depth))
    vafs, dp, reads = vafs[remove_zeros], dp[remove_zeros], reads[remove_zeros]

    # Generate parameter vector and outputs similar to fully stochastic positive selection simulation data
    parameters = [0.0, 0.0, 0.0, Nfinal, 0.0, 0.0, nclonal, depth, noise, rho, 0.0, 0.0, 0, 0.0, "none", alt_reads]
    n_detectable = 0
    beneficial, deleterious, neutral, total, subclonal = mutationMetrics(vafs, ["N" for i in 1:length(vafs)], [0.0 for i in 1:length(vafs)], f_min = f_min, f_max = f_max)

    return n_detectable, [vafs, reads, dp], [0.0, 0.0, 0.0, 0.0, 0.0], [beneficial, deleterious, neutral, total, subclonal], parameters
end
