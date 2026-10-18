---
title: Irreducibility testing over finite fields
summary: Rabin's test, Ben-Or's test and the distinct-degree test for irreducibility over $\mathbb{F}_q$, their costs on irreducible and on random inputs, how FLINT combines them, and the construction of irreducible polynomials.
complexity: For $f$ of degree $n$ over $\mathbb{F}_q$: Ben-Or and the distinct-degree test $O(\sqrt{n}\, C(n) + n M(n) + M(n) \log q)$ in the worst case, much less on random inputs; Rabin's test with a Frobenius ladder $O(\log n \cdot (1 + \omega(n))\, C(n) + M(n) \log q)$.
order: 78
---

Testing whether $f \in \mathbb{F}_q[x]$ of degree $n$ is irreducible is
cheaper than factoring it, and is the inner loop of the construction of
irreducible polynomials, which define the extension fields
$\mathbb{F}_{q^n}$. All tests are based on
[](algorithms/polynomials/factoring-finite-fields#thm-xqd): $x^{q^d} - x$
is the product of the monic irreducible polynomials whose degree divides
$d$. As in the factorisation algorithms, the powers $h_d = x^{q^d} \bmod
f$ are computed from $h_1 = x^q \bmod f$ either by exponentiation by $q$
or by modular composition, $h_{i + j} = h_i(h_j) \bmod f$
([](algorithms/polynomials/factoring-finite-fields#lem-frob)). $M(n)$
and $C(n)$ denote the costs of multiplication and of modular composition.

## Rabin's test {#rabin}

::: theorem Rabin [@Rab1980b] {#thm-rabin}
Let $f \in \mathbb{F}_q[x]$ have degree $n \ge 1$. Then $f$ is irreducible
if and only if
$$
x^{q^n} \equiv x \pmod f \quad \text{and} \quad \gcd\bigl(f, x^{q^{n/r}} -
x\bigr) = 1 \text{ for every prime } r \mid n.
$$
:::

::: proof
If $f$ is irreducible, it divides $x^{q^n} - x$, and it divides
$x^{q^{n/r}} - x$ for no $r$ because $n \nmid n/r$. Conversely, the first
condition says that $f$ divides the squarefree polynomial $x^{q^n} - x$,
so $f$ is squarefree and all its irreducible factors have degrees
dividing $n$. If $f$ had a factor of degree $d < n$, $d \mid n$, then $d$
would divide $n/r$ for some prime $r \mid n$, and that factor would
divide $\gcd(f, x^{q^{n/r}} - x)$.
:::

The test needs $h_n$ and $h_{n/r}$ for the $\omega(n)$ prime divisors $r$
of $n$. `gr_poly_is_irreducible_rabin` computes the ladder $L_i =
x^{q^{2^i}} \bmod f = L_{i-1}(L_{i-1}) \bmod f$ for $2^i \le n$, and then
each $h_e$ as a composition of the $L_i$ along the binary expansion of
$e$. This takes $O(\log n)$ compositions for the ladder and for each of
the $1 + \omega(n)$ powers, independently of $q$ apart from the first
step $L_0 = h_1$, which costs $O(M(n) \log q)$. The test does the same
amount of work whether or not $f$ is irreducible.

`nmod_poly_is_irreducible_rabin` is an older implementation of the same
test that computes each $h_e$ by $e$ exponentiations by $q$, at a cost of
$O(n M(n) \log q)$.

## Ben-Or's test and the distinct-degree test {#ben-or}

::: theorem Ben-Or [@BenOr1981] {#thm-ben-or}
$f$ of degree $n \ge 1$ is irreducible if and only if $\gcd(f, x^{q^i} -
x) = 1$ for $1 \le i \le n/2$.
:::

::: proof
A reducible $f$ has an irreducible factor of degree $d \le n/2$, which
divides $x^{q^d} - x$. Conversely, an irreducible $f$ of degree $n$
divides $x^{q^i} - x$ only if $n \mid i$.
:::

Ben-Or's test computes $h_1, h_2, \ldots$ one after the other and stops at
the first nontrivial gcd. Its worst case, reached for irreducible $f$, is
$n/2$ steps; but a random polynomial has a factor of small degree with
high probability: by an analogue of Mertens' theorem, the probability
that a random monic polynomial of degree $n$ has no irreducible factor of
degree at most $i$ tends to about $e^{-\gamma}/i$ for large $n$ and $i$, so
the expected number of steps on random input is small [@GaoPan1997]. Since a random polynomial of degree $n$
is irreducible with probability about $1/n$, this early exit is what
matters when searching for irreducible polynomials, and Gao and Panario
[@GaoPan1997] found Ben-Or's test much faster than Rabin's in that
setting. `gr_poly_is_irreducible_ben_or` implements it, with each step a
composition with the precomputed matrix of powers of $h_1$ when $q$ is
large, and an exponentiation by $q$ otherwise.

The *distinct-degree test* `gr_poly_is_irreducible_ddf` runs the coarse
phase of the baby-step giant-step distinct-degree factorisation
([](algorithms/polynomials/factoring-finite-fields#alg-ddf-bsgs)): it
computes $\ell \approx \sqrt{n/2}$ baby steps and then, for each interval
$(\ell j, \ell(j+1)]$ of degrees, one giant step and the gcd of $f$ with
the interval polynomial, stopping at the first nontrivial gcd. It is
Ben-Or's test with the degrees grouped into intervals, so it costs
$O(\sqrt{n}\, C(n) + n M(n))$ for irreducible inputs and also exits early
for most reducible ones, though only after the $\ell$ baby steps.

Both tests first check squarefreeness with one gcd ($\gcd(f, f') = 1$),
which is cheap and rejects some inputs immediately.

## FLINT's strategy {#strategy}

`gr_poly_is_irreducible`, which `nmod_poly_is_irreducible`,
`fmpz_mod_poly_is_irreducible` and the `fq_*_poly_is_irreducible`
functions wrap, proceeds as follows for $f$ of degree $n \ge 2$:

1. Over a prime field $\mathbb{F}_p$ with $p \le \max(200, 2(n + 1))$,
   evaluate $f$ at all $p$ elements and return false if there is a root.
   For small $p$ this rejects a large fraction of random inputs at a cost
   of $O(np)$.
2. If $n < 600$, return the result of the distinct-degree test.
3. Otherwise check squarefreeness, then compute $\gcd(f, x^{q^i} - x)$
   for $i \le 4$ (with exponentiations) and return false if one is
   nontrivial, and finally return the result of Rabin's test.

For large irreducible inputs Rabin's test with the Frobenius ladder needs
$O(\log n)$ compositions against $O(\sqrt{n})$ for the distinct-degree
test, so the crossover is at a moderate degree; the filter in step 3
keeps the early exit for random inputs, whose smallest factor usually
has small degree. The generic functions return `T_UNKNOWN` when they
cannot decide (for example over a context that is not known to be a
finite field); the type-specific wrappers throw an exception in that
case, and, unlike the generic functions, consider constants
irreducible.

## Constructing irreducible polynomials {#construction}

An irreducible polynomial of degree $n$ over $\mathbb{F}_q$ is found by
testing random monic polynomials: there are about $q^n / n$ of them among
the $q^n$ monic polynomials, so about $n$ trials are needed on average
(`nmod_poly_randtest_irreducible`, `fmpz_mod_poly_randtest_irreducible`).
For finite field arithmetic sparse polynomials are preferable, because
reduction modulo a trinomial or pentanomial is cheap;
`nmod_poly_minimal_irreducible` searches binomials, trinomials,
tetranomials and pentanomials in a fixed order with
`nmod_poly_is_irreducible_ddf`. Conway polynomials, which make the
embeddings between different extension fields compatible, cannot be
computed efficiently in general and are taken from a table
([@Lüb2004]; `_nmod_poly_conway`). The finite field contexts such as
`fq_nmod_ctx_init_ui` use the Conway polynomial when the table contains
it and the sparse polynomial of `nmod_poly_minimal_irreducible`
otherwise. The choice of defining polynomials for the finite field
types is discussed in [](algorithms/finite-fields/defining-polynomials).

::: functions gr_poly_is_irreducible gr_poly_is_irreducible_ddf gr_poly_is_irreducible_rabin gr_poly_is_irreducible_ben_or nmod_poly_is_irreducible nmod_poly_is_irreducible_ddf nmod_poly_is_irreducible_rabin fmpz_mod_poly_is_irreducible fmpz_mod_poly_is_irreducible_ddf fmpz_mod_poly_is_irreducible_rabin :::
