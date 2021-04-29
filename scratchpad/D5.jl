
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using VinbergsAlgorithmNF

using Hecke

using LinearAlgebra

using CoxeterDiagrams

using ToggleableAsserts
toggle(true)

VA = VinbergsAlgorithmNF

⊕ = VA.Lat.:(⊕)   
bd = VA.Lat.block_diag

A = VA.Lat.A
B = VA.Lat.B
D = VA.Lat.D
E = VA.Lat.E
II = VA.Lat.I

ℚ,a = Hecke.rationals_as_number_field()

mat = -1 .* II(1) ⊕ D(5)
vd = VinbergData(ℚ,mat)

roots_at_distance_zero = VA.roots_at_distance_zero(vd)
@assert roots_at_distance_zero == [
    [0, 1, 0, 0, 0, 0],
    [0, -1, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, -1, 0, 0, 0],
    [0, -1, -1, 0, 0, 0],
    [0, 1, 1, 1, 0, 0],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, -1, -1, 0, 0],
    [0, -1, -1, -1, 0, 0],
    [0, 0, 0, -1, 0, 0],
    [0, 1, 1, 1, 1, 0],
    [0, 0, 1, 1, 1, 0],
    [0, 0, 0, 1, 1, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, -1, -1, -1, 0],
    [0, -1, -1, -1, -1, 0],
    [0, 0, 0, -1, -1, 0],
    [0, 0, 0, 0, -1, 0],
    [0, 1, 1, 2, 1, 1],
    [0, 0, 1, 2, 1, 1],
    [0, 1, 2, 2, 1, 1],
    [0, 1, 1, 1, 1, 1],
    [0, 0, 1, 1, 1, 1],
    [0, 0, 0, 1, 1, 1],
    [0, 1, 1, 1, 0, 1],
    [0, 0, 1, 1, 0, 1],
    [0, 0, 0, 1, 0, 1],
    [0, 0, 0, 0, 0, 1],
    [0, 0, -1, -1, -1, -1],
    [0, -1, -1, -1, -1, -1],
    [0, 0, 0, -1, -1, -1],
    [0, 0, -1, -2, -1, -1],
    [0, -1, -1, -2, -1, -1],
    [0, -1, -2, -2, -1, -1],
    [0, 0, -1, -1, 0, -1],
    [0, -1, -1, -1, 0, -1],
    [0, 0, 0, -1, 0, -1],
    [0, 0, 0, 0, 0, -1],
    [0, 2, 2, 2, 1, 1],
    [0, 0, 2, 2, 1, 1],
    [0, 0, 0, 2, 1, 1],
    [0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, -1, 1],
    [0, 0, -2, -2, -1, -1],
    [0, -2, -2, -2, -1, -1],
    [0, 0, 0, -2, -1, -1],
    [0, 0, 0, 0, -1, -1],
    [0, 0, 0, 0, 1, -1]
]

cone_roots = VA.cone_roots(vd)
@assert cone_roots == [
    [0, -1, 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0],
    [0, 0, 0, -1, 0, 0],
    [0, 0, 0, 0, -1, 1],
    [0, 0, 0, 0, 0, -1]
]

expected_cone_roots = [
    [0, -1, 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0],
    [0, 0, 0, -1, 0, 0],
    [0, 0, 0, 0, -1, 0],
    [0, 0, 0, 0,  0,-1]
]

for expected_root in expected_cone_roots
    @assert !VA.is_necessary_halfspace(vd.gram_matrix.entries,cone_roots,expected_root) || !VA.is_necessary_halfspace(vd.gram_matrix.entries,cone_roots,-expected_root)
end
