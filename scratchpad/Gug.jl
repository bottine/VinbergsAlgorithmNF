using Hecke
using VinbergsAlgorithmNF

VA = VinbergsAlgorithmNF

# Define our field
Qx,x = Hecke.QQ["x"]
f = 8*x^3 + 4*x^2 - 4*x - 1     # minimal poly of cos(2Ο/7)
K,a = Hecke.NumberField(f, "a")

# simplify its internal representation or whatever (otherwise Hecke has some troubles doing computations)
L, mL = simplify(K)

# ring of algebraic integers πͺ_L
OL = Hecke.maximal_order(L)

Οβ = inv(mL)(a)  # = cos(6Ο/7)
Οβ = 2Οβ^2-1     # = cos(2Ο/7)
Οβ = 2Οβ^2-1     # = cos(4Ο/7)


# check that we really have the values we want
@assert VA.approx(Οβ) β cos(2Ο/7) 
@assert VA.approx(Οβ) β cos(4Ο/7)
@assert VA.approx(Οβ) β cos(6Ο/7)

# In Bugaenko, it's [Ο,Οβ,Οβ] (so, we're compatible with his naming)
# To match Guglielmetti we have aΟβ+bΟβ+cΟβ =  [a//2,b//2,c//2] (because you have to multiply by 2 to get elements of πͺ_L

# Helper function, to pass from Guglielmetti's representation to an element of L
gug_to_phi(vβ,vβ,vβ) = 2Οβ*vβ + 2Οβ*vβ + 2Οβ*vβ

# The β€-module rep of an element in πͺ_LΒ can be obtained by `.coordinates` or `coordinate(β¦)` which I think is better
# We build the matrix allowing to go from Guglielmetti's rep to Hecke's internal one.
OL_change_of_basis = Rational{Int}.(Int.(hcat([OL(2Οβ).coordinates,OL(2Οβ).coordinates,OL(2Οβ).coordinates]...)))

T = OL_change_of_basis
S = inv(OL_change_of_basis)

# We Check that it's the right matrix 
@assert S*(OL(2Οβ).coordinates) == [1,0,0]
@assert S*(OL(2Οβ).coordinates) == [0,1,0]
@assert S*(OL(2Οβ).coordinates) == [0,0,1]

# We define the diagonal of the form we want to look at.
# This is using Guglielmetti's representation
form = [
    gug_to_phi(-1,0,0),
    1,
    gug_to_phi(-2,-5,-7)
]

# We check that it's the right one using another representation (better safe than sorry)
@assert form == [-2Οβ,1,-4Οβ-10Οβ-14Οβ]

# The roots given by AlVin
gug_roots = [
    [0,-1,0],
    [0,0,-1],
    [1,gug_to_phi(0,0,-1),0],
    [gug_to_phi(0,-1,-2),gug_to_phi(0,0,-1),gug_to_phi(1,0,0)],
    [gug_to_phi(-1,-2,-4),0,gug_to_phi(0,-1,-1)],
]

# And another way to write them
@assert gug_roots == [[0,-1,0],[0,0,-1],[1,-2Οβ,0],[-2Οβ-4Οβ,-2Οβ,2Οβ],[-2Οβ-4Οβ-8Οβ,0,-2Οβ-2Οβ]]

# The roots my implementation finds, modulo some unit changes
my_roots = [[0, -1, 0],[0, 0, 2*(2Οβ)^2 + (2*Οβ) - 5],[2*(2Οβ)^2 + 2Οβ - 4, (2Οβ)^2 - 2, 0],[-2*(2Οβ)^2 - 4*Οβ + 5, 0, -3*(2Οβ)^2 - 4*Οβ + 7]]

# My roots and AlVin's one agree up to the third one (colinear takes care of rescalings):
@assert all(VA.colinear(my,gug) for (my,gug) in zip(my_roots[1:3],gug_roots[1:3]))

# Written in the "internal" representation
my_roots_Hecke_basis = [(coeff -> (coordinates(OL(coeff)))).(root) for root in my_roots]
# Written w.r.t. the basis used by Guglielmetti
my_roots_gug_basis = [[S * coeff for coeff in root] for root in my_roots_Hecke_basis] 

# Some sanity checks
for (root,Hecke_root,gug_root) in zip(my_roots,my_roots_Hecke_basis,my_roots_gug_basis)
    for (coeff,Hecke_coeff,gug_coeff) in zip(root,Hecke_root,gug_root)
        @assert Hecke_coeff == T*gug_coeff && S*Hecke_coeff == gug_coeff
        @assert gug_to_phi(gug_coeff...) == coeff
    end
end


println("My Roots in Guglielmetti's representation:")
display(my_roots_gug_basis)
println("------------------------------------------")



# Set up the algo's base data
vd = VA.VinbergData(L,VA.diagm(L,form))

# Check that it chooses the right basepoint
@assert VA.basepoint(vd) == [1,0,0]

# Run the algo
(status,(roots,dict,diagram)) = VA.next_n_roots!(vd,n=10)
@assert status == true

# To view the diagram, uncomment the following lines (and change the destination path to whatever you want):
# io = open("./nice_picture_of_a_coxeter_diagram.png","w")
# show(io, MIME"image/png"(), diagram)
# close(io)


# Two helper functions to check that the algo's output is reasonably what is expected (the first computes the distance to basepoint and orders them; the second computes the gram coefficients and orders them)
inc_dists(vd,roots) = sort(map(r -> VA.fake_dist_to_basepoint(vd,r),roots))
gram_coeffs(vd,roots) = sort([VA.Gram_coeff(vd.quad_space,rβ,rβ) for rβ in roots for rβ in roots])

# We check that the distances and gram coeffs computed are as expected (here we should compare the Coxeter diagram, but I haven't implemented any diagram isomorphism checking -- this is a "reasonable" substitute)
@assert inc_dists(vd,my_roots) == inc_dists(vd,roots)
@assert gram_coeffs(vd,my_roots) == gram_coeffs(vd,roots)
#@assert all(VA.colinear(my_r,r) for (my_r,r) in zip(my_roots,roots)) # does not work because the algo sometimes orders things differentlyβ¦ this should be remedied somehow

println("Distances to basepoint (GUG):")
for root in gug_roots
    println(VA.approx(VA.fake_dist_to_basepoint(vd,root)))
end
println("and inner products:")
println([VA.approx(VA.times(vd,rβ,rβ)) for rβ in gug_roots for rβ in gug_roots if rββ rβ])

println()
println("Distances to basepoint (MY):")
for root in my_roots
    println(VA.approx(VA.fake_dist_to_basepoint(vd,root)))
end

println("and inner products:")
println([VA.approx(VA.times(vd,rβ,rβ)) for rβ in my_roots for rβ in my_roots if rββ rβ])


