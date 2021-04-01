#=

@testset "AffineConstraints" begin
    K,a = Hecke.rationals_as_number_field()
    @test VA.interval_for_k_j(VA.no_constraints(),1) == (nothing,nothing)
    @test VA.interval_for_k_j(VA.no_constraints(),2) == (nothing,nothing)
    
    bot_left = ([K.([1,0]),K.([0,1])],K.([0,0]),[1,2])
    # 1x+0y ≤ 0
    # 0x+1y ≤ 0
    @test VA.interval_for_k_j(bot_left,1) == (nothing,0)
    @test VA.interval_for_k_j(bot_left,2) == (nothing,0)
    
    bot_right = ([K.([1,0]),K.([0,-1])],K.([0,0]),[1,2])
    @test VA.interval_for_k_j(bot_right,1) == (nothing,0)
    @test VA.interval_for_k_j(bot_right,2) == (0,nothing)
    
    above_x_eq_y = ([K.([1,-1])],K.([0]),[2])
    @test VA.interval_for_k_j(above_x_eq_y,1) == (nothing,nothing)
    @test VA.interval_for_k_j(above_x_eq_y,2) == (0,nothing)
    
    constraint = ([K.([1,-1]),K.([-1,-1])],K.([0,0]),[2,2])
    @test VA.interval_for_k_j(constraint,1) == (nothing,nothing) 
    @test VA.interval_for_k_j(constraint,2) == (0,nothing)
    
    constraint = ([K.([1,1]),K.([-1,-1])],K.([0,0]),[2,2]) # this is just a line
    @test VA.interval_for_k_j(constraint,1) == (nothing,nothing) # can't deduce a constraint on x
    @test VA.interval_for_k_j(constraint,2) == (0,0)       # but what remains constrains y
    
    constraint = ([K.([1,1]),K.([-1,-1])],K.([-1,-1]),[2,2]) # this is empty 
    @test VA.interval_for_k_j(constraint,1) == (nothing,nothing) # can't deduce a constraint on x!! because it's not "obvious" 
    @test VA.interval_for_k_j(constraint,2) == (1,-1)       # but what remains constrains y

end
=#
