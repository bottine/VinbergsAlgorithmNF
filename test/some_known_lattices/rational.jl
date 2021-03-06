rational_lattices = [
    (
        "Sasha_1a",
        [2 -1 0 0; -1 2 -1 0; 0 -1 2 0; 0 0 0 -1],
        [[-1,0,1,0],[0,-1,0,0],[0,0,-1,0],[1,1,1,1]]
    ),
    (
        "Sasha_1b",
        [2 -1 0 0; -1 2 -1 0; 0 -1 2 0; 0 0 0 -2],
        [[-1,0,1,0],[0,-1,0,0],[0,0,-1,0],[1,2,1,1],[3,2,1,2]]
    ),
    (
        "Sasha_1c",
        [2 -1 0 0; -1 2 -1 0; 0 -1 2 0; 0 0 0 -10], 
        [[-1,0,1,0],[0,-1,0,0],[0,0,-1,0],[5,5,5,2],[3,2,1,1],[5,10,5,3],[5,6,3,2]]
    ),
    (
        "Sasha_3",
        [-5 0 0 0; 0 1 0 0; 0 0 2 -1; 0 0 -1 2],
        [[0,-1,0,0],[0,0,-1,1],[0,0,0,-1],[1,1,2,1],[2,5,0,0],[1,2,1,1],[3,0,5,5]]
    ),
    (
        "Vin07_L19",
        [-8 0 0 0; 0 2 1 0; 0 1 2 1; 0 0 1 2],
        [[0, -1, 2, -1],[0, 0, -1, 1],[0, 0, 0, -1],[1, 2, 0, 2],[1, 3, -2, 1]]
    ),
    (
        "Vin07_L4",
        [0 3 0 0; 3 0 0 0; 0 0 2 0; 0 0 0 2],
        [[-1, -1, 0, 0], [0, 0, -1, 1], [0, 0, 0, -1], [1, 0, 1, 0], [2, -2, 3, 3]]
    ),
    (
        "Vin07_L6",
        [0 1 0 0; 1 0 0 0; 0 0 2 0; 0 0 0 8],
        [[-1, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1], [1, 0, 1, 0], [4, 0, 0, 1], [2, -2, 1, 1]],
    ),
    (
        "Vin07_L17",
        [0 1 0 0; 1 0 0 0; 0 0 2 1; 0 0 1 14],
        [[-1, -1, 0, 0], [0, 0, -1, 0], [0, 0, -1, 2], [1, 0, 1, 0], [3, -2, 1, -1], [27, 0, 1, -2], [6, -1, 1, -1]]
    ),
    (
        "BP18_12",
        [0 1 0 0; 1 0 0 0; 0 0 36 0; 0 0 0 6],
        [[-1, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1], [0, 3, 0, 1], [0, 18, 1, 0], [-3, 6, 1, 1], [-4, 4, 1, 0], [-2, 8, 1, 0], [-8, 8, 1, 4], [-2, 10, 1, 1], [-12, 12, 2, 5], [-36, 36, 7, 12], [-16, 22, 3, 8], [-54, 72, 11, 24], [-14, 22, 3, 7]]
    ),
    (
        "BP_D4D4_1",
        [2 0 1 0 0 0 0 0 0; 0 2 -1 0 0 0 0 0 0; 1 -1 2 -1 0 0 0 0 0; 0 0 -1 2 0 0 0 0 0; 0 0 0 0 2 0 1 0 0; 0 0 0 0 0 2 -1 0 0; 0 0 0 0 1 -1 2 -1 0; 0 0 0 0 0 0 -1 2 0; 0 0 0 0 0 0 0 0 -1],
        [[-1, 2, 2, 1, 0, 0, 0, 0, 0], [0, -1, 0, 1, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 2, 2, 1, 0], [0, 0, 0, 0, 0, -1, 0, 1, 0], [0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 1, 0, 0, 0, 1], [1, 0, 0, 0, 0, 0, 0, 0, 1], [2, -1, -2, -1, 2, -1, -2, -1, 2]]
    ),
    (
        "diag_-1,2,6,6",
        LinearAlgebra.diagm([-1,2,6,6]),
        [[0, -1, 0, 0], [0, 0, -1, 1], [0, 0, 0, -1], [1, 1, 0, 0], [2, 1, 1, 0], [3, 0, 1, 1]]
    ),
    (
        "diag_-2,1,1,1,1,1,1,1,1,1,1,1,1",
        LinearAlgebra.diagm([-2,1,1,1,1,1,1,1,1,1,1,1,1]),
        [[0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1], [1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0], [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
    ),
    (
        "diag_-20,2,2,3",
        LinearAlgebra.diagm([-20,2,2,3]),
        [[0, -1, 1, 0], [0, 0, -1, 0], [0, 0, 0, -1], [2, 5, 5, 0], [3, 10, 0, 0], [1, 2, 1, 2], [3, 9, 3, 2], [3, 6, 6, 4], [3, 0, 0, 8], [1, 3, 0, 1], [6, 15, 0, 10], [3, 8, 2, 4], [9, 15, 0, 20], [3, 3, 3, 7], [4, 10, 5, 5], [8, 5, 5, 20], [5, 6, 0, 12], [9, 9, 3, 22], [6, 5, 0, 15], [3, 4, 1, 7], [17, 20, 10, 40]]
    )

]


@testset "Some reflective lattices over ℚ" begin
    K,a = Hecke.rationals_as_number_field()
    
    for (name,matrix,known_roots) in rational_lattices
        @testset "$name" begin
            compare_vinberg_outputs(K,matrix,known_roots)
        end
    end

end

