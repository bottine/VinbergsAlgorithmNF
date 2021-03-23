Qx,x = Hecke.QQ["x"]
f = 8*x^3 + 4*x^2 - 4*x - 1
K,a = Hecke.NumberField(f, "a") 

# a = cos(6π/7)
# b = cos(2π/7)
# c = cos(4π/7

b = 2*a^2-1 # a is not equal to cos(2π/7) but b is
c = 2*b^2-1


lattices = [
    ( 
        "diag_-2ϕ,1,1,1",
        VA.diagm(K,[-2b,1,1,1]),
        [[0, -1, 1, 0], [0, 0, -1, 1], [0, 0, 0, -1], [4*a^2, 4*a^2 - 4*a - 1, 0, 0], [12*a^2 - 2*a - 1, 8*a^2 - 2*a - 1, 8*a^2 - 2*a - 1, 8*a^2 - 2*a - 1]]
    ),
    (
        "diag_-2ϕ,1,-", # Found in Guglielmetti page 132
        VA.diagm(K,[-2b,1,-4b-10c-14a]),
        ["Here should be 5 vectors"]

    )
]


@testset "Some reflective lattices over ℚ(cos(2π/7))" begin
    
    for (name,matrix,known_roots) in lattices
        @testset "$name" begin
            compare_vinberg_outputs(K,matrix,known_roots)
        end
    end

end

