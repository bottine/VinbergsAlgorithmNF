# Correctness

* Make cone root computation exact

# Efficiency

* Reduce allocations needed in `update_constraints` and `_extend_root_stem`? 
* Cone roots: One can fix an initial halfspace and only consider subsequent halfspaces when they have acute angle with the first one → should divide the time to find cone roots.
* Check which of the conditions on stem extension are most useful and which are costly and order them accordingly.
* Whenever possible, use the diagonal form to compute inner products
* Profile and optimize whatever can be.

# Quality

* Unit tests everywhere
* More examples to check good behaviour
* More asserts
* More non-diagonal lattices over bad fields
* Implement Coxeter diagram isomorphism test

# Theory

* Given roots at distance zero, can one algebraically find a set of cone roots?
* Can one construct a set of cone roots without even being given the roots at distance zero?
* If not, can one at least *check* that a set of cone roots is minimal?
* Find proof of the condition that root lengths divide 2*last_invariant_factor?

# Organisation

* Clean up `Manifest.toml`
* Disassemble `some_lattice.jl` and reorganize into different file(s) (probably move everything to `quad_form.jl`)
  Moreover, add lattice constructors for all irreducible spherical root systems, and add way to mix and match (like `lattice_playground.jl` in the old project)
* Stop using `Hecke.QuadraticForm`s since we only actually use `inner_product()` which is probably easier to do by hand.
* Allow running "short" and "long" tests independently (like they do in `Hecke.jl`)
* In the interval management part of `_extend_root_stem`, we have to take care of `α_j` everywhere because the constraints are given by the diagonalized versions of the previous roots and the constraints are given by the diagonal form: We should move the `α_j`s to the constraints to get something cleaner in `_extend_root_stem`
