using Pkg
Pkg.activate(".")

using VinbergsAlgorithmNF
using Hecke
using LinearAlgebra
using Logging
using CoxeterDiagrams
using ToggleableAsserts
toggle(false)

VA = VinbergsAlgorithmNF


⊕ = VA.Lat.:(⊕)    

K,a = Hecke.quadratic_field(2) # a is -√2

# The following lattices LX correspond (up to typos) to the lattices L(X) in arXiv:2003.11944v2

L3 = [
    2   -1+a;
    -1+a   2
] ⊕ VA.Lat.I(2) # should be reflective (Prop 7.2)
 
L4 = [
    2   -1+2a;
    -1+2a   2
] ⊕ VA.Lat.I(2) # should be reflective (Prop 7.3)
 
L9 = [
    2   -1  0 ;
    -1   2  -1 ;
    0   -1   a
] ⊕ VA.Lat.I(1) # should _not_ be reflective (Prop 7.4)
 
L10 = [
    2   -1+a;
    -1+a   2
]  ⊕  (2-a) .* VA.Lat.I(1) ⊕ VA.Lat.I(1) # should _not_ (Prop 7.5)
 
L13 = [
     2   0   0  -1;
     0   2  -1  -1;
     0  -1   2   a;
    -1  -1   a   2
] # should _not_ (Prop 7.6)
 
L14 = [
     2   0   0    -1;
     0   2  -1    -1;
     0  -1   2     a-1;
    -1  -1   a-1   2
] # should _not_ (Prop 7.6)

L15 = [
     2   0  -1    -1;
     0   2  -1     a;
    -1  -1   2     a-1;
    -1   a   a-1   2
] # should _not_ (Prop 7.6)



io = open("log.txt", "w+")
logger = SimpleLogger(io)
global_logger(logger)


# Let's see what we get:
for lat in [L3,L4,L9,L10,L13,L14,L15]
    println()
    println("Lattice       : ")
    display(lat)
    println() 
    vd = VinbergData(K,lat)
    println("Basepoint     : ", VA.basepoint(vd))
    out = @timed VA.next_n_roots!(vd,n=10)
    (status,(roots,dict,diagram)) = out.value
    time = out.time

    println("time          : ", time)
    println("finite volume : ", status)
    println("roots         : ")
    display(roots)
    println()
end

close(io)


