
K,a = Hecke.quadratic_field(5)
ϕ = (-a+1)//2

sqrt_5_lattices = [
    (
        "diag_-ϕ,1,1,1,2",
        VA.diagm(K,[-ϕ,1,1,1,2]),
        [[0, -1, 1, 0, 0], [0, 0, -1, 1, 0], [0, 0, 0, -1, 0], [0, 0, 0, 0, -1], [1//2*a + 3//2, 1, 0, 0, 0], [1, -1//2*a - 1//2, 0, 0, 1], [-1//2*a + 1//2, 1, 1, 1, 1]]

    ),
    (
        "diag_-ϕ,1,1,1,1",
        VA.diagm(K,[-ϕ,1,1,1,1]),
        [[0, -1, 1, 0, 0], [0, 0, -1, 1, 0], [0, 0, 0, -1, 1], [0, 0, 0, 0, -1], [1//2*a + 3//2, 1, 0, 0, 0]]
    ),
    (
        "diag_-ϕ,1,1,1,1,1",
        VA.diagm(K,[-ϕ,1,1,1,1,1]),
        [[0, -1, 1, 0, 0, 0], [0, 0, -1, 1, 0, 0], [0, 0, 0, -1, 1, 0], [0, 0, 0, 0, -1, 1], [0, 0, 0, 0, 0, -1], [1//2*a + 3//2, 1, 0, 0, 0, 0], [-1//2*a + 1//2, 1, 1, 1, 1, 1]]
    ),
    #( # Takes too long, haven't even let it finish yet!
    #    "diag_-ϕ,1,1,1,1,1,1",
    #    LinearAlgebra.diagm([-ϕ,1,1,1,1,1,1]),
    #),
    (
        "Belolipetsky15",
        [ 2 -1  0  0; 
         -1  2 -ϕ  0; 
          0 -ϕ  2 -1;
          0  0 -1  2],
        [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [-1//2*a + 3//2, -a + 3, -3//2*a + 3//2, -1//2*a - 1//2]]

    )
]

@testset "Some reflective lattices over ℚ(√5)" begin
    
    for (name,matrix,known_roots) in sqrt_5_lattices
        @testset "$name" begin
            compare_vinberg_outputs(K,matrix,known_roots)
        end
    end

end

