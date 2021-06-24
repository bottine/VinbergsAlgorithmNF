
@testset "Infinite order symmetries for simple diagonal forms over ℚ" begin
    K,a = Hecke.rationals_as_number_field()


    # cf Guglielmetti's thesis p. 143
    cases = [
        ([11,-1,1],true),
        ([11,-1,1,1],false),
        ([13,-1,1],true),
        ([13,-1,1,1],false),
        ([-7,1,1,1],true),
        ([-7,1,1,1,1],false),
        ([-13,1,1],true),
        ([-13,1,1,1],false),
        ([-14,1,1],true),
        ([-14,1,1,1],false),
        ([-17,1,1,1],true),
        ([-17,1,1,1,1],false),
    ]

    for (diag,reflective) in cases
        mat = VA.Lat.block_diag(diag...)
        vd = VinbergData(K,mat)
        (st,(root,dict,das)) = VA.next_n_roots!(vd,n=10)
        @test st == reflective
        @test VA.inf_ord_sym2(vd,root,das) != reflective
        @test VA.inf_ord_sym3(vd,root,das) != reflective
    end

end


@testset "Infinite order symmetries for forms over ℚ(√2)" begin
    K,a = Hecke.quadratic_field(2)


    # cf Bogachev's paper 
    cases = [
        (
            [
                 2  0  0  -1;
                 0  2 -1  -1;
                 0 -1  2   a-1;
                -1 -1  a-1 2
            ],
            false
        ),
    ]

    for (mat,reflective) in cases
        vd = VinbergData(K,mat)
        (st,(root,dict,das)) = VA.next_n_roots!(vd,n=10)
        @test st == reflective
        @test VA.inf_ord_sym2(vd,root,das) != reflective
    end

end


