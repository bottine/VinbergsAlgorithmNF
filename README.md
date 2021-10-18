# 𝐕𝐢𝐧𝐧𝐲 𝐭𝐡𝐞 𝐜𝐨𝐱𝐛𝐨𝐲

WIP implementation of Vinberg's algorithm over number fields (with PID algebraic integers).
In `Julia`, using [`Hecke.jl`](https://github.com/thofma/Hecke.jl) for all number theoretic computations.
Based on/inspired by

* [Guglielmetti's `AlVin`](https://github.com/rgugliel/AlVin)
* [Bogachev & Perepechko's `VinAl`](https://github.com/aperep/vinal)

The performance is reasonable on small enough examples: Finding Bugaenko's 8 dimensional polyhedron (from the lattice -(√5+1) ⟂ E₈) took <10mn last time I tried.
The program chokes on some seemingly easy examples (say the matrix obtained by replacing 0 by -(√5+1)/2 at positions (6,7) and (7,6) in the matrix E₆ ⟂ E₆).
This is due to to the way the program handles non-diagonal matrices.

### *Learning to herd Coxeter diagrams in the Hyperbolic Plains since 2021*
### *Vinny is now working the fields (slowly)*

### Howto

Assuming you have julia installed (A recent version is needed, e.g. it works 1.6.3 but not 1.5.3), it should go like this:

First, clone the repo:

``` 
git clone https://github.com/bottine/VinbergsAlgorithmNF
```

Move to the folder:

```
cd VinbergsAlgorithmNF
```

Launch julia:

```
julia
```

From julia, activate/instantiate the environment (whatever that means):

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

This should install whatever is necessary.
To run the tests (which takes some time and is not necessary):

```julia
Pkg.test()
```

Otherwise, to run an example:

```julia
using Hecke
using VinbergsAlgorithmNF; VA = VinbergsAlgorithmNF
```

And now one can try to run the algo (staying in the julia REPL):
E.g, to run 

#### the small example from Belolipetsky survey:

```julia
K,a = Hecke.quadratic_field(5) # a is -√5
ϕ = (a-1)//2
Bel = [2 -1 0 0; -1 2  -ϕ 0; 0 -ϕ 2 -1; 0 0 -1 2]
vd = VA.VinbergData(K,Bel)
(status, (roots,dict,das)) = VA.next_n_roots!(vd,n=10)
```

Now `status` tells whether the roots found form a finite volume diagram, `roots` is the roots found and `dict,das` are internal structures used for the iteration.

**Normally**, calling `show(das)` in jupyter or the like should show the Coxeter diagram defined by the roots.

To continue running the algorithm where we stopped (which is not necessary here, since there are less than 10 roots), call

```julia
(status, (roots,dict,das) = VA.next_n_roots!(vd,roots,dict,das,n=10)
```

#### On ℚ

To use ℚ as the field of definition, we need to let

```julia
K,a = Hecke.rationals_as_number_field()
```

#### Remarks

* The code is still messy, buggy and slow.
* The best way to understand how the functions above work is probably to just look at the code for now.
* Questions welcome.
* See the file(s?) in `scratchpad/` for a (some?) example(s?).




