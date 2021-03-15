using Test
using VinbergsAlgorithmNF
using Hecke

VA = VinbergsAlgorithmNF

@testset "Known cases" begin

    @testset "Rationals" begin
        K,a = Hecke.rationals_as_number_field()

        (status,output) = VA.all_in_one(K.([-1,2,6,6]),10)
        @test status == true
        @test all(VA.colinear(v₁,v₂) for (v₁,v₂) in zip(output[1],[[0, -1, 0, 0], [0, 0, -1, 1], [0, 0, 0, -1], [1, 1, 0, 0], [2, 1, 1, 0], [3, 0, 1, 1]]))
	

        (status,output) = VA.all_in_one(K.([-2,1,1,1,1]),10)
        @test status == true
        @test all(VA.colinear(v₁,v₂) for (v₁,v₂) in zip(output[1],[[0, -1, 1, 0, 0], [0, 0, -1, 1, 0], [0, 0, 0, -1, 1], [0, 0, 0, 0, -1], [1, 2, 0, 0, 0], [1, 1, 1, 1, 1]]))
        
    end

    @testset "ℚ(√5)" begin
        K,a = quadratic_field(5)
        (status,output) = VA.all_in_one(K.([(-1+a)//2,1,1,1,2]),10)
        @test status == true
        @test all(VA.colinear(v₁,v₂) for (v₁,v₂) in zip(output[1],[[0, -1, 1, 0, 0], [0, 0, -1, 1, 0], [0, 0, 0, -1, 0], [0, 0, 0, 0, -1], [1//2*a + 3//2, 1, 0, 0, 0], [1, -1//2*a - 1//2, 0, 0, 1], [-1//2*a + 1//2, 1, 1, 1, 1]]))  
    end

end
