
K,a = Hecke.quadratic_field(5)
ϕ = (-a+1)//2
⊕ = VA.Lat.:(⊕)

sqrt_5_lattices = [
    (
        "Bugaenko 7dim first example (https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/poly/7-all.pdf)",
        -2ϕ .* VA.Lat.I(1) ⊕ VA.Lat.E(7),
        [[0, -1, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, -1], [1, -a + 1, -a + 1, -3//2*a + 3//2, -2*a + 2, -3//2*a + 3//2, -a + 1, -1//2*a + 1//2], [-1//2*a + 1//2, -a + 1, -3//2*a + 3//2, -2*a + 2, -3*a + 3, -5//2*a + 5//2, -2*a + 2, -a + 1], [-1//2*a + 3//2, -a + 3, -3//2*a + 7//2, -2*a + 5, -3*a + 7, -5//2*a + 11//2, -2*a + 4, -3//2*a + 5//2], [-a + 2, -2*a + 4, -7//2*a + 13//2, -4*a + 8, -6*a + 12, -9//2*a + 19//2, -3*a + 7, -3//2*a + 9//2]]    
    ),
    (
        "Bugaenko 8dim first example (https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/poly/8.pdf)",
        -2ϕ .* VA.Lat.I(1) ⊕ VA.Lat.E(8),
        [[0, -1, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1], [-1//2*a - 1//2, 2, 3, 4, 6, 5, 4, 3, 2], [1, 4, 5, 7, 10, 8, 6, 4, 2], [-a + 2, -3*a + 7, -9//2*a + 23//2, -6*a + 14, -9*a + 21, -15//2*a + 33//2, -6*a + 12, -4*a + 8, -2*a + 4]] 
    ),   
]

@testset "Bugaenko's big reflective lattices over ℚ(√5)" begin
    
    for (name,matrix,known_roots) in sqrt_5_lattices
        @testset "$name" begin
            compare_vinberg_outputs(K,matrix,known_roots)
        end
    end

end

