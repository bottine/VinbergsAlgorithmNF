
K,a = Hecke.quadratic_field(2)
sqrt_2_lattices = [
    (
        "diag_-ϕ,1,1,1,2",
        VA.diagm(K,[a-1,1,1,1,1]),
        [[0, -1, 1, 0, 0], [0, 0, -1, 1, 0], [0, 0, 0, -1, 1], [0, 0, 0, 0, -1], [1, -a + 1, 0, 0, 0], [-a + 1, -a + 1, -a + 1, -a + 1, 0], [-a + 2, -a + 2, -a + 1, -a + 1, -a + 1]]
    ),
    
]

@testset "Some reflective lattices over ℚ(√2)" begin
    
    for (name,matrix,known_roots) in sqrt_2_lattices
        @testset "$name" begin
            compare_vinberg_outputs(K,matrix,known_roots)
        end
    end

end

