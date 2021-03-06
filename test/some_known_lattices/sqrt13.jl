
K,a = Hecke.quadratic_field(13)
ϕ = (-a+1)//2

lattices = [
    (
        "SimilartoBelolipetsky15",
        [2 -1 0 0; -1 2  -ϕ 0; 0 -ϕ 2 -1; 0 0 -1 2],
        [[-1, 0, 0, 0], [0, -1, 0, 0], [1, 2, -1//2*a - 1//2, -2], [-1, -2, 1//2*a + 1//2, 1], [0, 0, -1, 0], [1, 1//2*a - 3//2, 1//2*a - 3//2, 1]]
    )
]

@testset "Some reflective lattices over ℚ(√13)" begin
    
    for (name,matrix,known_roots) in lattices
        @testset "$name" begin
            compare_vinberg_outputs(K,matrix,known_roots)
        end
    end

end

