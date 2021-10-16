# Roots

The goal being to enumerate the mirror hyperplanes of hyperbolic reflections contained in ``Γ_r``.

***

A reflection ``ρ`` of ``V`` has the form

```math
    x ↦ x - 2\frac{xr}{r²}r
```

where ``r∈V`` is the *root* of the reflection.
The reflection mirror ``H_ρ`` is the set of fixed points of ``ρ``, or equivalently the set ``r^⟂`` of elements ``x`` of ``V`` with ``xr = 0``.

We are only interested in reflections along roots whose mirrors intersect ``ℍ^n``, or equivalently, roots ``r`` satisfying ``r²<0``.
Furthermore, for the reflection ``ρ = ρ_r`` (along root ``r``) to lie in ``𝒪^L``, ``r`` has to have entries in ``K`` (i.e. all coordinates of ``r`` must lie in our field, otherwise there is no chance to preserve the lattice).
Since ``R`` is principal, we can rescale ``r`` and assume that all its coordinate lie in ``R`` itself.
In that case, let's say that ``r`` is *integral*.
By principality again, we can also assume that the gcd of the entries of ``r`` is a unit in ``R``, in which case we say that ``r`` is *primitive*.
Finally, let us say that ``r`` satisfies the crystallographic condition if ``ρ_r∈𝒪^L``.

Fix ``r`` integral satisfying the crystallographic condition, i.e. with

```math
    ∀x∈L:\ x - 2\frac{xr}{r²}r ∈ L.
```

This condition is clearly equivalent to

```math
    2\frac{e_ir}{r²}r ∈ L.
```

for all basis vector ``e_i`` of ``L``.
Assume now that ``r`` is primitive.
Since its entries are coprime, the condition

```math
    2\frac{e_ir}{r²}r ∈ L
```

simply means that ``r²`` divides ``2e_ir``.
Thus, a *primitive* vector ``r`` satisfies the crystallographic condition iff ``r² | 2e_ir`` for all ``e_i``.

From now on, we call a root only those vectors that are integral, primitive and satisfy the crystallographic condition.
Call also ``r²`` the *length* of ``r``.

## Root lengths

An important fact due to Vinberg:

**Lemma.** 
If ``r`` is a root, then ``r²`` divides twice the last invariant factor of ``Q``.

*Proof.*
Let 

```math
  L^* ≔ \{x ∈ V : q(x,y) ∈ R\ ∀  y ∈ L\}
```

Since ``Q`` is integral, ``L⊆L^*`` and the decomposition of ``L^*/L`` into cyclic submodules (which exists since ``R`` is principal) is exactly the decomposition given by the Smith Normal Form of ``Q``.
In particular, if ``f`` is the last invariant factor of ``Q``, then ``fx ∈ L`` for all ``x∈L^*``.
The crystallographic condition for ``r`` states exactly that ``\frac{2r}{r²}∈L^*``.
It follows that ``\frac{2f}{r²}r ∈ L``, and using primitivity again, that ``r²`` divides ``2f``. 

The conclusion is that, forgetting units, the number of possible 

## Units

By fancy number theory results, in our case the group of units ``U`` of ``R`` is a finitely generated (abelian) group.
Thus, ``U/U²`` is a finite group, which by the way Hecke can deal with.

If ``r`` is a root and ``u`` is a unit, then ``ur`` is *still* a root with exact same induced reflection.
The lengths of ``r`` and ``ur`` differ by ``u²``.
Let ``V`` be a set of representatives in ``U`` of ``U/U²``.
Say that two roots ``r,s`` are *equivalent* if ``s = ur`` for some unit ``u``.

**Fact.**
Up to sign, any equivalence class of roots is represented by a unique root ``r`` with ``r²∈UD``, where ``D`` is a set of representative divisors of ``2f`` (up to units) and ``f`` is the last invariant factor of ``Q``.

*Proof.*
Pretty much by construction (see my thesis, p.70).

**Note.**
Obviously, if ``r`` is a root, so is ``-r``, of exact same length (hence the “up to sign” in the fact above).
The difference between ``r`` and ``-r`` is that they correspond to the two possible orientations of the reflection mirror ``H_r``, so the lack of unicity here is not problematic.

In conclusion: Enumerating *oriented* hyperplanes of hyperbolic reflection mirrors of ``Γ_r`` amounts to enumerating roots of length contained in `UD`, and those cover all possible hyperplanes **canonically**.

## Root lengths part 2

A possible length ``l∈UD`` may turn up to have non-positive Galois conjugates.

If ``r`` is a root, and ``σ`` any non-trivial Galois embedding, then ``(σr)²`` must be strictly positive, since the non-trivial conjugates of ``Q`` are positive-definite.
It follows that if ``σ(l)`` is non-positive, it cannot be the length of a root, and we may discard it.

## Roots by distance to a basepoint.

If ``r`` is a root, and ``p₀∈V`` with ``p₀²<0`` represents a choice of basepoint in ``ℍ^n``, then the distance between ``H_r`` and ``p₀`` in ``ℍ^n`` is given by:

```math
  \sinh⁻¹ \sqrt{-\frac{(p₀r)²}{p₀²r²}}.
```

Moreover, ``p₀`` lies in the negative halfspace defined by ``r`` iff ``p₀r≤0`` (the case ``p₀r=0`` means that ``p₀`` is contained in the reflection mirror of ``r``).

Assuming that ``p₀`` is fixed and *does* lie in the negative halfspace of ``r``, the distance (in hyperbolic space) between ``p₀`` and ``H_r`` increases monotonically with

```math
  \frac{p₀r}{r²}. 
```
