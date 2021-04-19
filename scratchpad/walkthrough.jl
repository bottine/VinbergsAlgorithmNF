# This script should be callable just with `julia this.jl` from the command line 
# Or simply pasted to julia's repl

using Pkg
# First, activate the project: this should load all dependencies, etc
Pkg.activate(".")
# `instantiate` actually ensures that the dependency are exactly at the same state than on my computer
Pkg.instantiate()

# Load the module
using VinbergsAlgorithmNF

# Load Hecke, which deals with everything number theoretic
using Hecke

# And LinearAlgebra
using LinearAlgebra

# And CoxeterDiagrams, which contains the algo to check co{compactness,finiteness}
using CoxeterDiagrams

# **Important**: if you want the code to run fast enough, always call `toggle(false)` first.
# I peppered the code with asserts at some performance critical places, and calling `toggle(false)` disables those asserts.
using ToggleableAsserts
toggle(false)

# This is just a shortcut to not have to write `VinbergsAlgorithmNF` in full everytime
VA = VinbergsAlgorithmNF

# This is "orthogonal sum"
⊕ = VA.Lat.:(⊕)   
# This is "orthogonal sum" with as many summand as you want
bd = VA.Lat.block_diag

# These are the usual lattices (defined in `src/some_lattices.jl`)
A = VA.Lat.A
B = VA.Lat.B
D = VA.Lat.D
E = VA.Lat.E
II = VA.Lat.I

# Define the field ℚ(√5).
# WARNING: a == -√5 (by default the generator is the negative one)
K,a = quadratic_field(5)


Bugaenko7 = (a-1) .* II(1) ⊕ E(7) # The notation `.*` in `(a-1) .* II(1)` means that we multiply each entry in the matrix `II(1)` by `a-1`
Bugaenko8 = (a-1) .* II(1) ⊕ E(8)
Belolipetsky4 = A(2) ⊕ A(2); Belolipetsky4[3,2] = (a-1)//2; Belolipetsky4[2,3] = (a-1)//2
some_other = bd((a-1)//2 .* II(1),A(4),D(4)) 

# WARNING: Due to my code for Coxeter diagrams being messy, please run only one instance of the algorithm at a time!
#          Note also that for now the algo only handles up to 256 vertices!


# Initialize the needed data: need to specify the matrix and the field
vd = VinbergData(K,Bugaenko7)

# `VA.next_n_roots!` looks for at least `n` roots _after_ the cone roots.
# `status` is a boolean telling whether we got a finite volume polyhedron
# `roots` is the list of roots we got
# `dict` contains data used in the algorithm (should not touch it)
# `das` contains data related to the Coxeter diagram
(status,(roots,dict,diagram)) = VA.next_n_roots!(vd,n=10)

# Here we got all the roots we need, so no need to go further.
println("Are we finite volume? $status")
display(roots)

vd = VinbergData(K,Bugaenko8)
(status,out) = VA.next_n_roots!(vd,n=1)

# For Bugaenko8, there are three roots after the cone roots
println("Are we finite volume? $status")
display(out[1])

# So, we need to call the algo again: we feed it the necessary "runtime" data:
(status,out) = VA.next_n_roots!(vd,out...,n=1)
println("Are we finite volume? $status")
display(out[1])

# The last root takes more time to find 
# (status,out) = VA.next_n_roots(vd,out...,n=1)

cone_roots_Bugaenko8 = [
    [0, -1, 0, 0, 0, 0, 0, 0, 0], 
    [0, 0, -1, 0, 0, 0, 0, 0, 0], 
    [0, 0, 0, -1, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, -1, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, -1, 0, 0, 0], 
    [0, 0, 0, 0, 0, 0, -1, 0, 0], 
    [0, 0, 0, 0, 0, 0, 0, -1, 0], 
    [0, 0, 0, 0, 0, 0, 0, 0, -1]
]
cone_roots_Bugaenko8 = [K.(r) for r in cone_roots_Bugaenko8] # need to force the coefficients to be in the field

# If we know the cone roots, we can also provide them beforehand:
(status,out) = VA.next_n_roots!(vd,cone_roots_Bugaenko8,n=1)
# WARNING: `next_n_roots!` adds the new roots to the one it was fed (in this case `cone_roots_Bugaenko8`), so that now `cone_roots_Bugaenko8` contains more roots than before

# WARNING: for now, I haven't implemented a way to feed the code an arbitrary sequence of roots, so if we want to stop the algo and resume it, we need to do it as above.
#          thus, when providing cone roots, one need to be sure that they really are all the cone roots: nothing more, nothing less!
