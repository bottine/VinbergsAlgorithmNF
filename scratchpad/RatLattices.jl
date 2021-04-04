using VinbergsAlgorithmNF
using Hecke
using LinearAlgebra
using Logging
using CoxeterDiagrams
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


# Let's see what we get:
for lat in vcat([S1,S2,B1,B2,B3,C1,C2,C3],[C4(n) for n in 3:7])
    println()
    println("Lattice       : ")
    display(lat)
    
    vd = VinbergData(ℚ,lat)
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

print("\n"^10)
print("COMPARING OUTPUTS!")
println()

# PR as in Proposed roots

# Both outputs of the sage VinAl version
PR_S1 =  [[0, 0, -1, 0], [0, 0, 1, -2], [1, 1, 0, 0], [-1, 0, 1, 0], [-2, 2, 0, 1], [-19, 0, -1, 2], [-4, 1, 0, 1]]
PR_S2 = [[0, 0, -1, 0], [0, 0, 1, -2], [1, -1, 0, 0], [0, 1, 1, 0], [2, 2, 0, 1], [0, 19, -1, 2], [1, 4, 0, 1]]

# Simpy output: but the first root appears again negatively!
PR_B1 = [[-3, -5, -1, 1], [-1, 0, 0, 0], [0, -1, -1, 0], [3, 5, 1, -1], [0, -1, 0, 0]]

# Sage output
PR_B2 = [[-1, -5, -1, 1], [1, 0, -1, 0], [0, 1, 1, 0], [1, 8, 2, -1], [0, 1, 0, 0]]

# Simpy output
PR_B3 = [[-1, -1, -1, 1], [-1, 0, 0, 0], [0, -1, 1, 0], [3, 1, 0, 0], [2, 2, 1, -1]]

for (lat,roots) in [(S1,PR_S1),(S2,PR_S2),(B1,PR_B1),(B2,PR_B2),(B3,PR_B3)]
    vd = VinbergData(ℚ,lat)
    
    (status,(my_roots,dict,diagram)) = VA.next_n_roots!(vd,n=20)
    @assert status # finishes with finite volume polytope

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


