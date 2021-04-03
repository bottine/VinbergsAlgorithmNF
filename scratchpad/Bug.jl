using VinbergsAlgorithmNF
using Hecke
using LinearAlgebra
VA = VinbergsAlgorithmNF

K,a = quadratic_field(5) # a is -√5

# The 9-dimensional lattice in Bugaenko's paper
Bug = VA.Lat.:(⊕)((a-1).*VA.Lat.I(1), K.(VA.Lat.E(8)))

# Takes long but should work
vd = VinbergData(K,Bug)
(status,(roots,dict,diagram)) = @time VA.next_n_roots!(vd,n=20)

@assert status == true
display(roots)
display(diagram)

