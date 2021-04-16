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


# Comparison with AlVin's SAGE and SimPy versions

⊕ = VA.Lat.:(⊕)    


ℚ,a = Hecke.rationals_as_number_field() 

U = VA.Lat.U()
D = VA.Lat.D(4)
N = (-1).*VA.Lat.I(1)

# S as in Strange
S1 = U ⊕ [2 1; 1 10]
S2 = (-1).*U ⊕ [2 1; 1 10]

# B as in Bad
B1 = [1 -1  1 -1;
     -1  1 -3 -1;
      1 -3  9 -3;
     -1 -1 -3  1]

B2 = [9 -3  3 -3;
     -3  1 -3 -1;
      3 -3  9 -3;
     -3 -1 -3  1]

B3 = [1 -1 -1 -1;
     -1  1 -1 -1;
     -1 -1  1 -1;
     -1 -1 -1  1]

# C as in Comparison
C1 = N ⊕ D 
C2 = N ⊕ D ⊕ D
C3 = N ⊕ D ⊕ D ⊕ D
C4(n) = (3 .* N) ⊕ VA.Lat.I(n)

io = open("log.txt", "w+")
logger = SimpleLogger(io)
global_logger(logger)


#=
# Let's see what we get:
for (lat,basepoint) in vcat([(S1,[-1,1,0,0]),(S2,nothing),(B2,nothing),(B3,nothing),(C1,nothing),(C2,nothing),(C3,nothing)],[(C4(n),nothing) for n in 3:7])
    println()
    println("Lattice       : ")
    display(lat)
    println() 
    vd = VinbergData(ℚ,lat,basepoint)
    println("Basepoint     : ", VA.basepoint(vd))
    out = @timed VA.next_n_roots!(vd,n=20)
    (status,(roots,dict,diagram)) = out.value
    time = out.time

    println("time          : ", time)
    println("finite volume : ", status)
    println("roots         : ")
    display(roots)
    println()
end
=#
print("COMPARING OUTPUTS:")
println()

# PR as in Proposed roots

# Both outputs of the sage VinAl version
PR_S1 =  [[0, 0, -1, 0], [0, 0, 1, -2], [1, 1, 0, 0], [-1, 0, 1, 0], [-2, 2, 0, 1], [-19, 0, -1, 2], [-4, 1, 0, 1]]
PR_S2 = [[0, 0, -1, 0], [0, 0, 1, -2], [1, -1, 0, 0], [0, 1, 1, 0], [2, 2, 0, 1], [0, 19, -1, 2], [1, 4, 0, 1]]


# Sage output
PR_B2 = [[-1, -5, -1, 1], [1, 0, -1, 0], [0, 1, 1, 0], [1, 8, 2, -1], [0, 1, 0, 0]]

# Simpy output
PR_B3 = [[-1, -1, -1, 1], [-1, 0, 0, 0], [0, -1, 1, 0], [3, 1, 0, 0], [2, 2, 1, -1]]

for (lat,roots,basepoint) in [(S1,PR_S1,[-1,1,0,0]),(S2,PR_S2,nothing),(B2,PR_B2,nothing),(B3,PR_B3,nothing)]
    println()
    println("Lattice       : ")
    display(lat)
    println() 

    # First run it on the cone roots contained in `roots` -> assuming same basepoint, the results should then be _exactly_ identical.
    vd = VinbergData(ℚ,lat,basepoint)
    proposed_cone_roots = deepcopy([ℚ.(r) for r in roots if VA.fake_dist_to_basepoint(vd,r)==0])
    out =  @timed VA.next_n_roots!(vd,deepcopy(proposed_cone_roots),n=20)
    (status,(my_roots,dict,diagram)) = out.value
    time = out.time
    @assert status # finishes with finite volume polytope
    @assert roots == my_roots

    # Then run without specifying roots -> the diagrams should still agree
    out =  @timed VA.next_n_roots!(vd,n=20)
    (status,(my_roots,dict,diagram)) = out.value
    time = out.time


    if !CoxeterDiagrams.is_isom(VA.Coxeter_matrix(vd,my_roots),VA.Coxeter_matrix(vd,roots))
        println("different Coxeter diagrams for:")
        display(lat)
        println()
        println("given roots:")
        for r in roots
            println("$r         at dist $(VA.fake_dist_to_basepoint(vd,r))")
        end
        println("my roots:")
        display(my_roots)
        for r in my_roots
            println("$r         at dist $(VA.fake_dist_to_basepoint(vd,r))")
        end
    end

end

close(io)


