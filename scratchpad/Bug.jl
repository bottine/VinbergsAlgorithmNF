using VinbergsAlgorithmNF
using Hecke
using LinearAlgebra


K,a = quadratic_field(5) # a is -√5

# The 9-dimensional lattice in Bugaenko's paper
Bug = VA.:(⊕)((a-1).*K.(VA.I_), K.(VA.gram_E8))

# Takes long but should work
vd = VinbergData(K,Bug)
(status,(roots,dict,diagram)) = VA.next_n_roots!(vd,n=20)

@assert status == true
display(roots)
display(diagram)

