# ππ’π§π§π² π­π‘π ππ¨π±ππ¨π²

WIP implementation of Vinberg's algorithm over number fields (with PIDΒ algebraic integers).
In `Julia`, using [`Hecke.jl`](https://github.com/thofma/Hecke.jl) for all number theoretic computations.
Based on/inspired by

* [Guglielmetti's `AlVin`](https://github.com/rgugliel/AlVin)
* [Bogachev & Perepechko's `VinAl`](https://github.com/aperep/vinal)

The performance is reasonable on small enough examples: Finding Bugaenko's 8 dimensional polyhedron (from the lattice -(β5+1) β Eβ) took <10mn last time I tried.
The program chokes on some seemingly easy examples (say the matrix obtained by replacing 0 by -(β5+1)/2 at positions (6,7) and (7,6) in the matrix Eβ β Eβ).
This is due to to the way the program handles non-diagonal matrices.

### *Learning to herd Coxeter diagrams in the Hyperbolic Plains since 2021*
### *Vinny is now working the fields (slowly)*

### Howto

Assuming you have julia installed, it should go like this:

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
using Pkg;Β Pkg.activate("."); Pkg.instantiate()
```

This should install whatever is necessary.
Then

```julia
using Hecke
using VinbergsAlgorithmNF; VA = VinbergsAlgorithmNF
```

And now one can try to run the algo (staying in the julia REPL):
E.g, to run 

#### the small example from Belolipetsky survey:

```julia
K,a = Hecke.quadratic_field(5) # a is -β5
Ο = (a-1)//2
Bel = [2 -1 0 0; -1 2Β  -Ο 0; 0 -Ο 2 -1; 0 0 -1 2]
vd = VA.VinbergData(K,Bel)
(status, (roots,dict,das)) = VA.next_n_roots!(vd,n=10)
```

Now `status` tells whether the roots found form a finite volume diagram, `roots` is the roots found and `dict,das` are internal structures used for the iteration.

**Normally**, calling `show(das)` in jupyter or the like should show the Coxeter diagram defined by the roots.

To continue running the algorithm where we stopped (which is not necessary here, since there are less than 10 roots), call

```julia
(status, (roots,dict,das) = VA.next_n_roots!(vd,roots,dict,das,n=10)
```

#### On β

To use β as the field of definition, we need to let

```julia
K,a = Hecke.rationals_as_number_field()
```

#### Remarks

* The code is still messy, buggy and slow.
* The best way to understand how the functions above work is probably to just look at the code for now.
* Questions welcome.
* See the file(s?) in `scratchpad/` for a (some?) example(s?).




