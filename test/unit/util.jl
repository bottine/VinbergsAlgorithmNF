

@testset "Diagonalize" begin
    
    ℚ,a = Hecke.rationals_as_number_field()
    ℤ = maximal_order(ℚ)
    
    k = 10

    for (d,r) in [(3,150),(4,100),(5,50),(6,20)], round in 1:r
           dim = d
           M = rand(-(k-d):(k-d),dim,dim)
           M = M + M'
           v = rand(-(k-d):(k-d),dim); while v'*M*v == 0 v = rand(-5:5,dim) end
           (D,P) = VA.diagonalize(ℤ,M,v)
           @test P'*M*P == D; 
           @test P[:,1] == v
    end
end

