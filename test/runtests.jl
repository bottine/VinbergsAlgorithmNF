using Test
using VinbergsAlgorithmNF
using Hecke
using LinearAlgebra
using CoxeterDiagrams
VA = VinbergsAlgorithmNF
using ToggleableAsserts
# Unit tests, eventually:




# Helper functions
inc_dists(vd,roots) = sort(map(r -> VA.fake_dist_to_basepoint(vd,r),roots))
gram_coeffs(vd,roots) = sort([VA.Gram_coeff(vd.gram_matrix,r₁,r₂) for r₁ in roots for r₂ in roots])

# test "similarity" of known_roots for reflective lattice with output of running VA.next_n_roots.
# This checks that the distances to the basepoints agree and that the list of Gram coefficients agree, but it might be that the Coxeter diagrams still are different -> this should be remedied.
function compare_vinberg_outputs(matrix,known_roots)
    vd = VA.VinbergData(matrix)
    (status,(roots_found,dict,das)) = VA.next_n_roots!(vd,n=length(known_roots))
    @test status == true
    @test length(roots_found) == length(known_roots)
    @test inc_dists(vd,roots_found) == inc_dists(vd,known_roots) 
    @test gram_coeffs(vd,roots_found) == gram_coeffs(vd,known_roots)
    @test CoxeterDiagrams.is_isom(VA.Coxeter_matrix(vd,roots_found),VA.Coxeter_matrix(vd,known_roots))
end


@testset "Known lattices" begin
    include("some_known_lattices/rational.jl")
    include("some_known_lattices/Bianchi.jl")
end
