using Test
using VinbergsAlgorithmNF
using Hecke
using LinearAlgebra

VA = VinbergsAlgorithmNF

# Helper functions
inc_dists(vd,roots) = sort(map(r -> VA.fake_dist_to_basepoint(vd,r),roots))
function compare_vinberg_outputs(K,matrix,known_roots)
    vd = VA.VinbergData(K,matrix)
    (status,(roots_found,dict,das)) = VA.next_n_roots!(vd,n=length(known_roots))
    @test status == true
    @test length(roots_found) == length(known_roots)
    @test inc_dists(vd,roots_found) == inc_dists(vd,known_roots)   
end


include("test_known_rational_lattices.jl")
include("test_known_sqrt2_lattices.jl")
include("test_known_sqrt5_lattices.jl")
include("test_known_sqrt13_lattices.jl")
# include("test_known_cos2piover7_lattices.jl") # strange errors with our choice of irred poly. New versions of Hecke should make it doable
