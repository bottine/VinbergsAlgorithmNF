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

        (status,output) = VA.all_in_one(K.([-2,1,1,1,1,1,1,1,1,1,1,1,1]),10)
        @test status = true
        @test all(VA.colinear(v₁,v₂) for (v₁,v₂) in zip(output[1],[[0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1], [1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0], [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]))
        
        (status,output) = VA.all_in_one(K.([-20,2,2,3]),20)
        @test status = true
        @test all(VA.colinear(v₁,v₂) for (v₁,v₂) in zip(output[1],[[0, -1, 1, 0], [0, 0, -1, 0], [0, 0, 0, -1], [2, 5, 5, 0], [3, 10, 0, 0], [1, 2, 1, 2], [3, 9, 3, 2], [3, 6, 6, 4], [3, 0, 0, 8], [1, 3, 0, 1], [6, 15, 0, 10], [3, 8, 2, 4], [9, 15, 0, 20], [3, 3, 3, 7], [4, 10, 5, 5], [8, 5, 5, 20], [5, 6, 0, 12], [9, 9, 3, 22], [6, 5, 0, 15], [3, 4, 1, 7], [17, 20, 10, 40]]))

        
    end

    @testset "ℚ(√5)" begin
        K,a = quadratic_field(5)
        (status,output) = VA.all_in_one(K.([(-1+a)//2,1,1,1,2]),10)
        @test status == true
        @test all(VA.colinear(v₁,v₂) for (v₁,v₂) in zip(output[1],[[0, -1, 1, 0, 0], [0, 0, -1, 1, 0], [0, 0, 0, -1, 0], [0, 0, 0, 0, -1], [1//2*a + 3//2, 1, 0, 0, 0], [1, -1//2*a - 1//2, 0, 0, 1], [-1//2*a + 1//2, 1, 1, 1, 1]]))  
    end

    @testset "ℚ(cos(2π/7))" begin
        Qx,x = Hecke.QQ["x"]
        f = x^3+(1//2)*x^2-(1//2)*x-1
        K,a = Hecke.NumberField(f, "a") 
        b = 2*a^2-1 # a is not equal to cos(2π/7) but b is
        (status,output) = VA.all_in_one(K.([-2b,1,1,1]),10)
        @test status == true
        @test all(VA.colinear(v₁,v₂) for (v₁,v₂) in zip(output[1],[[0, -1, 1, 0], [0, 0, -1, 1], [0, 0, 0, -1], [4*a^2, 4*a^2 - 4*a - 1, 0, 0], [12*a^2 - 2*a - 1, 8*a^2 - 2*a - 1, 8*a^2 - 2*a - 1, 8*a^2 - 2*a - 1]]))

        # Anything longer takes too long!
    end

end
