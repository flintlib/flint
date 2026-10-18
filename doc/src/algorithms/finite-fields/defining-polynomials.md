---
title: Defining polynomials
summary: Which irreducible polynomial defines $\mathbb{F}_{p^d}$ in FLINT: Conway polynomials from Lübeck's tables when available, otherwise polynomials with as few terms as possible (`fq_nmod`), random sparse irreducible polynomials (`fq`) or random primitive polynomials (`fq_zech`); with the facts about their existence that these choices rely on.
complexity: A table lookup for Conway polynomials; an expected $O(d)$ irreducibility tests for random search; primitivity testing requires the factorisation of $p^d - 1$.
order: 50
---

A field $\mathbb{F}_{p^d}$ is constructed as $\mathbb{F}_p[x]/(f)$ for an
irreducible $f$ of degree $d$. Any such $f$ gives an isomorphic field, but
the choice matters in practice:

- a *sparse* $f$, with few nonzero terms, makes the reduction after each
  multiplication cheap ([](algorithms/finite-fields/reduction));
- a *primitive* $f$, whose roots generate $\mathbb{F}_q^\times$, is
  required by the Zech logarithm representation and gives a known
  generator of the multiplicative group (used, for example, as the
  non-square in `fq_sqrt`);
- a *standard* $f$ makes results reproducible across systems and allows
  compatible embeddings between fields of different degrees.

## Counting irreducible and primitive polynomials {#counting}

The number of monic irreducible polynomials of degree $d$ over
$\mathbb{F}_p$ is $\frac{1}{d} \sum_{e \mid d} \mu(d/e) p^e$, which lies
between $p^d/d - 2p^{d/2}/d$ and $p^d/d$ [@LN1997]; so a random monic
polynomial of degree $d$ is irreducible with probability about $1/d$,
and random search with an irreducibility test
([](algorithms/polynomials/irreducibility)) needs about $d$
attempts. Of the $p^d - 1$ nonzero elements, $\varphi(p^d - 1)$ are
primitive, and each primitive polynomial has $d$ primitive roots, so
there are $\varphi(p^d - 1)/d$ primitive polynomials, roughly a fraction
$\varphi(p^d - 1)/(p^d - 1)$ of the irreducible ones, which is at least
$c/\log\log(p^d)$ for a constant $c > 0$.
Testing primitivity requires the factorisation of $p^d - 1$
(`fq_multiplicative_order`).

## Conway polynomials {#conway}

::: definition {#def-conway}
The Conway polynomial $C_{p,d}$ is the least monic primitive polynomial
$f$ of degree $d$ over $\mathbb{F}_p$, in a fixed ordering of
polynomials, such that for every proper divisor $m$ of $d$, the
polynomial $C_{p,m}(x^{(p^d - 1)/(p^m - 1)})$ is divisible by $f$.
:::

The ordering is lexicographic on the sequence
$(-f_{d-1}, f_{d-2}, -f_{d-3}, \ldots, (-1)^d f_0)$ of coefficients with
alternating signs, each reduced to an integer in $[0, p)$. The compatibility condition
says that if $\alpha_d$ is a root of $C_{p,d}$ and $m \mid d$, then
$\alpha_d^{(p^d - 1)/(p^m - 1)}$ is a root of $C_{p,m}$; so the roots
of the Conway polynomials of all degrees form a compatible system, and
the field $\mathbb{F}_p[x]/(C_{p,m})$ embeds into
$\mathbb{F}_p[x]/(C_{p,d})$ by $x \mapsto x^{(p^d-1)/(p^m-1)}$. Conway
polynomials were defined by Parker for the computation of modular
character tables and named after J. H. Conway; they are used by GAP,
Magma, Sage and other systems as the standard defining polynomials.
Computing them is expensive, since the definition requires a search
through primitive polynomials (with the factorisation of $p^d - 1$) and
compatibility tests; Heath and Loehr give algorithms [@HeaLoe2004], and
Lübeck maintains the reference tables [@Lüb2004].

FLINT stores a compressed copy of Lübeck's table in `nmod_poly`
(`_nmod_poly_conway`), for all primes $p < 109988$ and the degrees for
which a Conway polynomial is known: degrees up to 409 for $p = 2$, up to
about 200 for $p < 20$, and fewer for larger primes (up to 12 for
$p < 300$, 9 for $p < 3371$, 6 for $p < 11000$ and 4 above), in each case
with gaps. Since Conway polynomials are primitive, they are always
suitable for the Zech representation, but they are not sparse in
general.

## Sparse and minimal-weight polynomials {#sparse}

When no Conway polynomial is available, the modules choose sparse
polynomials, which minimise the cost of reduction.

**Binomials.** By a classical theorem [@LN1997], $x^d - a$ is irreducible
over $\mathbb{F}_p$ if and only if every prime divisor of $d$ divides the
order $e$ of $a$ but not $(p - 1)/e$, and $p \equiv 1 \pmod 4$ if
$4 \mid d$. In particular an irreducible binomial of degree $d$ exists
only if every prime factor of $d$ divides $p - 1$.

**Trinomials** $x^d + a x^k + b$ exist for most but not all $(p, d)$: by
Swan's theorem there is no irreducible trinomial of degree
$d \equiv 0 \pmod 8$ over $\mathbb{F}_2$ [@Swa1962], and other
obstructions are known. It is conjectured that an irreducible
tetranomial ($p > 2$) or pentanomial ($p = 2$) always exists.

`nmod_poly_minimal_irreducible` returns a monic irreducible polynomial
with the least possible number of terms, trying binomials, trinomials,
tetranomials and pentanomials in this order and returning the first one
in a fixed ordering, after sieving out candidates with small factors; it
is deterministic. `fmpz_mod_poly_randtest_sparse_irreducible` draws
random irreducible trinomials, then pentanomials, then random sparse
polynomials, for multiprecision $p$.

## Choices in FLINT {#choices}

| Module | Default defining polynomial |
|---|---|
| `fq` (`fq_ctx_init`) | Conway polynomial if $p < 109988$ and it is in the table; otherwise a random sparse irreducible polynomial (trinomial or pentanomial if one is found), with a fixed random state |
| `fq_nmod` (`fq_nmod_ctx_init_ui`) | Conway polynomial if in the table; otherwise the minimal-weight polynomial of `nmod_poly_minimal_irreducible` |
| `fq_zech` (`fq_zech_ctx_init_ui`) | Conway polynomial if in the table; otherwise a random monic primitive polynomial (`nmod_poly_randtest_monic_primitive`, which draws random irreducible polynomials until the root is primitive) |
| `fq_default` | that of the representation it chooses |

In all modules a polynomial can also be given explicitly
(`fq_ctx_init_modulus` and its analogues); it must be irreducible, and
primitive for `fq_zech`. The choices for `fq` and `fq_zech` without a
Conway polynomial are reproducible for a given FLINT version, but are
not standardised; programs that exchange elements with other software
should give the defining polynomial explicitly.

## Use in FLINT {#use}

- `_nmod_poly_conway` (table lookup) is used by
  `_fq_nmod_ctx_init_conway_ui`, `_fq_ctx_init_conway_ui` and
  `_fq_zech_ctx_init_conway_ui`, and through them by all default
  context constructors.
- `fq_nmod_ctx_init_minimal_weight_ui` uses
  `nmod_poly_minimal_irreducible`; `fq_ctx_init` uses
  `fmpz_mod_poly_randtest_sparse_irreducible`; `fq_zech_ctx_init_random_ui`
  uses `nmod_poly_randtest_monic_primitive`.
- `fq_poly_randtest_irreducible` constructs random irreducible
  polynomials over $\mathbb{F}_q$ itself by random search with a
  distinct-degree test.
