using CanEvolve
using Test

@testset "CanEvolve.jl" begin
	p, n = autoSimulation(log(2), 0.0, Nfinal = 1000, noise = "betabinomial", rho = 0.002, nsubclone_min = 1, nsubclone_max = 2, lower_cutoff = 0.1, upper_cutoff = 0.4)
	@test sum(p[3][1]) > 0.0
	@test sum(n[3][1]) == 0.0
end
