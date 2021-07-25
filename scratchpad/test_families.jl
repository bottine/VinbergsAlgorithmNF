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


sanitize(mat) = begin
   str = string(mat)
   str = replace(str, "Any[" => s"")
   str = replace(str, "]" => s"")
   str = replace(str, r"//" => s"_over_")
   return str
end



function run(io,K,mat,num_roots,mat_short=Nothing)
    vd = VinbergData(K,mat)
    (st,(roots,dict,das)) = VA.next_n_roots!(vd,n=2) 
    println(io) 
    println(io, "Field: $K")
    println(io, "Matrix: ($(mat_short===nothing ? "" : mat_short)):")
    show(io, MIME"text/plain"(), mat)
    println(io) 
    while !st && length(roots) < num_roots
        inf,info = VA.inf_ord_sym_n_plus_one_walls(vd,roots,das,:vertices)
        if inf
            println(io, "Finite volume:")
            println(io, !inf)
            println(io, "First roots:")
            show(io, MIME"text/plain"(), roots)
            println(io)
            println(io, "Root permutation:")
            println(io, info[1])
            println(io, "Isometry:")
            show(io, MIME"text/plain"(), info[2])
            println(io)
            return
        end 

        (st,(roots,dict,das)) = VA.next_n_roots!(vd,roots,dict,das,n=2) 
         
    end
   
    if st
        println(io, "Finite volume:")
        println(io, st)
        println(io, "Roots ($(length(roots))):")
        show(io, MIME"text/plain"(), roots)
        println(io)
    else
        @warn "INCONCLUSIVE!!!"
        println(io, "Inconclusive")
    end
end


Qx, x = PolynomialRing(FlintQQ, "x");


K5,a5 = Hecke.NumberField(x^2-5,"(-√5)")
ϕ5 = (a5-1)//2
rows_Bugaenko = [
        (K5,2ϕ5⊕E(6),20,"2ϕ⊕E₆"),
        (K5,ϕ5⊕2E(6),20,"ϕ⊕2E₆"),
        (K5,2ϕ5⊕E(7),20,"2ϕ⊕E₇"),
        (K5,ϕ5⊕2E(7),20,"ϕ⊕2E₇"),
        (K5,2ϕ5⊕E(8),20,"2ϕ⊕E₈"),
        (K5,ϕ5⊕2E(8),20,"ϕ⊕2E₈"),
]

rows_K5bis = [
    (K5,2ϕ5⊕D(4),20,"2ϕ⊕D₄"),
    (K5,2ϕ5⊕D(5),20,"2ϕ⊕D₅"),
    (K5,2ϕ5⊕D(6),20,"2ϕ⊕D₆"),
    (K5,2ϕ5⊕D(7),20,"2ϕ⊕D₇"),
    (K5,2ϕ5⊕A(4),20,"2ϕ⊕A₄"),
    (K5,2ϕ5⊕A(5),20,"2ϕ⊕A₅"),
    (K5,2ϕ5⊕A(6),20,"2ϕ⊕A₆"),
    (K5,2ϕ5⊕A(7),20,"2ϕ⊕A₇"),
    (K5,2ϕ5⊕B(4),20,"2ϕ⊕B₄"),
    (K5,2ϕ5⊕B(5),20,"2ϕ⊕B₅"),
    (K5,2ϕ5⊕B(6),20,"2ϕ⊕B₆"),
    (K5,2ϕ5⊕B(7),20,"2ϕ⊕B₇"),
]

ℚ,aℚ = Hecke.rationals_as_number_field()
rows_BP2 = [(ℚ,ℚ.(-k⊕A(3)),20,"-$k⊕A₃") for k in 1:18]
rows_BP1 = [(ℚ,ℚ.(-2⊕A(2)⊕II(k)),16+k,"-2⊕A₂⊕I($k)") for k in 0:9]

families = [
    ("moreK5", rows_K5bis),
    ("BP1",rows_BP1),
    ("BP2",rows_BP2),
    ("Bugaenko",rows_Bugaenko),
]

for (name, rows) in families
    open("./families/$name", "w") do io
        println("OPENED")
        println("======")
        println()
        for row in rows
            (K,mat,num_roots,mat_short) = row
            println()
            println("mat: ", mat_short)
            println()
            println(io)
            run(io,K,mat,num_roots,mat_short)
            println(io)
            flush(io)
        end
    end
end


