---
title: fmpz_mod_poly_factor.h
header: fmpz_mod_poly_factor.h
description: Factorisation, irreducibility and squarefreeness testing, and root finding for polynomials over $\Z/n\Z$ with multiword modulus.
category: Univariate polynomials
order: 9
---

This module factors polynomials over the prime field $\mathbb{F}_p$
represented as `fmpz_mod_poly_t` (for primes of any size), tests them for
irreducibility and squarefreeness, and finds their roots, over
$\mathbb{F}_p$ and also over $\Z/n\Z$ for composite $n$ with known
factorisation. It mirrors `nmod_poly_factor`, and as there, the
factorisation, irreducibility and root-finding functions are wrappers
around the generic finite-field code of `gr_poly`
(`gr_poly_factor_finite_field`, `gr_poly_is_irreducible`,
`gr_poly_roots_finite_field`). For moduli below $2^{\bits}$,
`nmod_poly_factor` is faster.

A complete factorisation consists of a squarefree factorisation and, for
each squarefree part, the Kaltofen–Shoup algorithm (baby-step giant-step
distinct-degree factorisation with modular composition, and equal-degree
splitting), except for polynomials of very small degree, where the
Cantor–Zassenhaus algorithm is used; see
[](nmod_poly_factor/factor) for the rule.

::: summary algorithms/polynomials/factoring-finite-fields :::

::: summary algorithms/polynomials/root-finding :::

## Types and conventions {#conventions}

::: type
fmpz_mod_poly_factor_struct
fmpz_mod_poly_factor_t

A polynomial in factored form. The structure has the fields `poly` (an
array of `alloc` polynomials, of which the first `num` are used), `exp`
(an array of `alloc` exponents of type `slong`), `num` and `alloc`.
Entry $i$ for $0 \le i < \mathit{num}$ stands for
$\mathtt{poly[i]}^{\mathtt{exp[i]}}$. An `fmpz_mod_poly_factor_t` is an
array of one `fmpz_mod_poly_factor_struct`, so it is passed by reference.
:::

- All functions take the modulus as an `fmpz_mod_ctx_t`, which must be
  the context of all polynomials involved.
- The factorisation, irreducibility, squarefreeness and root-finding
  functions require the modulus $p$ to be prime, except
  `fmpz_mod_poly_roots_factored`, `fmpz_mod_poly_is_squarefree_f` and
  `fmpz_mod_poly_is_irreducible_rabin_f`. This is not checked: a composite
  modulus gives wrong results or an exception.
- The factors computed are monic. Unlike `nmod_poly_factor`, the
  factorisation functions do not return the leading coefficient.
- The factors are stored in no particular order. The algorithms are
  randomised with a fixed seed, so the output is reproducible.
- Most functions overwrite the output factorisation, but
  `fmpz_mod_poly_factor_distinct_deg` and
  `fmpz_mod_poly_factor_equal_deg` append to it.
- Distinct-degree and equal-degree factorisation of polynomials of degree
  above about 200 use several threads when they have been made available
  with `flint_set_num_threads`.

## Memory management {#memory}

::: function
void fmpz_mod_poly_factor_init(fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)

Initialises `fac` as the empty factorisation, without allocating memory.
:::

::: function
void fmpz_mod_poly_factor_clear(fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)

Frees the memory used by `fac`, including all allocated polynomials.
:::

::: function
void fmpz_mod_poly_factor_realloc(fmpz_mod_poly_factor_t fac, slong alloc, const fmpz_mod_ctx_t ctx)

Changes the space of `fac` to exactly `alloc` factors. For `alloc` = 0
the factorisation is cleared and reinitialised as empty. If `alloc` is
smaller than the number of factors, the factors beyond the first `alloc`
are discarded.
:::

::: function
void fmpz_mod_poly_factor_fit_length(fmpz_mod_poly_factor_t fac, slong len, const fmpz_mod_ctx_t ctx)

Ensures that `fac` has space for at least `len` factors. If it has to
grow, it grows to at least twice its previous allocation.
:::

::: function
void fmpz_mod_poly_factor_swap(fmpz_mod_poly_factor_t a, fmpz_mod_poly_factor_t b, const fmpz_mod_ctx_t ctx)

Swaps `a` and `b` efficiently.
:::

## Basic manipulation {#manipulation}

::: function
void fmpz_mod_poly_factor_set(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)

Sets `res` to a copy of `fac`.
:::

::: function
void fmpz_mod_poly_factor_get_poly(fmpz_mod_poly_t a, const fmpz_mod_poly_factor_t b, slong i, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_factor_get_fmpz_mod_poly(fmpz_mod_poly_t z, fmpz_mod_poly_factor_t fac, slong i, const fmpz_mod_ctx_t ctx)

Set the output to the factor of index $i$, $0 \le i < \mathit{num}$. The
exponent is `exp[i]`. `fmpz_mod_poly_factor_get_fmpz_mod_poly` is an
older name, with a non-`const` factorisation argument.
:::

::: function
void fmpz_mod_poly_factor_insert(fmpz_mod_poly_factor_t fac, const fmpz_mod_poly_t poly, slong exp, const fmpz_mod_ctx_t ctx)

Multiplies the factorisation by $\mathit{poly}^{\mathit{exp}}$: if
`poly` is equal to one of the factors, `exp` is added to its exponent,
otherwise a copy of `poly` is appended with exponent `exp`. Constant
polynomials are ignored. The search is linear in the number of factors.
:::

::: function
void fmpz_mod_poly_factor_concat(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)

Multiplies `res` by the factorisation `fac`, by inserting each factor of
`fac` with `fmpz_mod_poly_factor_insert`. `res` and `fac` must be
distinct objects.
:::

::: function
void fmpz_mod_poly_factor_pow(fmpz_mod_poly_factor_t fac, slong exp, const fmpz_mod_ctx_t ctx)

Multiplies every exponent of `fac` by `exp`, which raises the factored
polynomial to the power `exp`.
:::

## Input and output {#io}

::: function
void fmpz_mod_poly_factor_print(const fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_factor_print_pretty(const fmpz_mod_poly_factor_t fac, const char * var, const fmpz_mod_ctx_t ctx)

Print the factors to standard output, one per line, as the polynomial
followed by ` ^ ` and the exponent, in the formats of
`fmpz_mod_poly_print` and of `fmpz_mod_poly_print_pretty` with the
variable name `var`.
:::

## Factorisation {#factorisation}

::: function
void fmpz_mod_poly_factor(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Sets `res` to the factorisation of `f` into monic irreducible
polynomials, with multiplicities. The leading coefficient of `f` is not
returned. A constant `f` gives the empty factorisation.

::: details
This is `gr_poly_factor_finite_field` over the `fmpz_mod` context:
deflation, squarefree factorisation, and then the Cantor–Zassenhaus
algorithm if the degree $n$ of the (deflated) input satisfies $n <
125/(\log_2 p)^{1.3}$, which for primes of 40 bits or more never
happens, and the Kaltofen–Shoup algorithm otherwise. See
[](algorithms/polynomials/factoring-finite-fields#choice).
:::
:::

The following example factors $x^{15} - 1$ modulo the Mersenne prime
$p = 2^{127} - 1$. Since $p \equiv 7 \pmod{15}$, the multiplicative order
of $p$ modulo 15 is 4, so the cyclotomic factor $\Phi_{15}$ of degree 8
splits into two factors of degree 4, $\Phi_5$ stays irreducible (the
order of $p$ modulo 5 is 4) and $\Phi_3$ splits into two linear factors,
which with $x - 1$ gives three roots. The factors are not sorted by
degree, so the example sorts the degrees before printing them.

```c example name=fmpz_mod_poly_factor_cyclotomic
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>

static int cmp(const void * a, const void * b)
{
    slong x = *(const slong *) a, y = *(const slong *) b;
    return (x > y) - (x < y);
}

int main(void)
{
    fmpz_t p;
    fmpz_mod_ctx_t ctx;
    fmpz_mod_poly_t f;
    fmpz_mod_poly_factor_t fac;
    slong i, deg[16];

    fmpz_init(p);
    fmpz_ui_pow_ui(p, 2, 127);
    fmpz_sub_ui(p, p, 1);                     /* p = 2^127 - 1 = 7 mod 15 */
    fmpz_mod_ctx_init(ctx, p);
    fmpz_mod_poly_init(f, ctx);
    fmpz_mod_poly_factor_init(fac, ctx);

    fmpz_mod_poly_set_coeff_ui(f, 15, 1, ctx);
    fmpz_mod_poly_set_coeff_si(f, 0, -1, ctx);      /* x^15 - 1 */

    fmpz_mod_poly_factor(fac, f, ctx);
    for (i = 0; i < fac->num; i++)
        deg[i] = fmpz_mod_poly_degree(fac->poly + i, ctx);
    qsort(deg, fac->num, sizeof(slong), cmp);
    flint_printf("degrees of the factors:");
    for (i = 0; i < fac->num; i++)
        flint_printf(" %wd", deg[i]);
    flint_printf("\n");

    fmpz_mod_poly_roots(fac, f, 0, ctx);
    flint_printf("number of roots: %wd\n", fac->num);

    fmpz_mod_poly_factor_clear(fac, ctx);
    fmpz_mod_poly_clear(f, ctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(p);
    return 0;
}
```
```output
degrees of the factors: 1 1 1 4 4 4
number of roots: 3
```

::: function
void fmpz_mod_poly_factor_cantor_zassenhaus(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_factor_berlekamp(fmpz_mod_poly_factor_t factors, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_factor_kaltofen_shoup(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Like `fmpz_mod_poly_factor`, but split the squarefree factors with the
Cantor–Zassenhaus algorithm, Berlekamp's algorithm or the Kaltofen–Shoup
algorithm respectively (deflation is still applied). The output is
overwritten.
:::

::: function
void fmpz_mod_poly_factor_squarefree(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Sets `res` to the squarefree factorisation of `f`: monic squarefree
polynomials $g_i$, pairwise coprime and with distinct exponents $e_i$,
such that $f = c \prod_i g_i^{e_i}$ with $c$ the leading coefficient of
$f$. A constant `f` gives the empty factorisation. This is
`gr_poly_factor_squarefree`; see [](algorithms/polynomials/squarefree).
:::

::: function
void fmpz_mod_poly_factor_distinct_deg(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t poly, slong * const * degs, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_factor_distinct_deg_threaded(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t poly, slong * const * degs, const fmpz_mod_ctx_t ctx)

Distinct-degree factorisation of a squarefree polynomial `poly` of
degree at least 1, which need not be monic. Appends to `res` monic
polynomials $f_1, \ldots, f_k$, each with exponent 1, such that $f_j$ is
the product of all monic irreducible factors of `poly` of one degree
$d_j$, and writes $d_j$ to `(*degs)[j - 1]`. The array `*degs` must have
room for $\deg(\mathit{poly})$ entries; the $d_j$ are distinct but not
sorted. If `poly` is not squarefree the output is meaningless.

::: details
This is `gr_poly_factor_distinct_deg`, the baby-step giant-step
algorithm of Kaltofen and Shoup [@KalSho1998]; see
[](algorithms/polynomials/factoring-finite-fields#ddf). It uses several
threads by itself when they are available, so the `_threaded` version is
an alias. A constant `poly` causes an exception.
:::
:::

::: function
void fmpz_mod_poly_factor_equal_deg(fmpz_mod_poly_factor_t factors, const fmpz_mod_poly_t pol, slong d, const fmpz_mod_ctx_t ctx)

Equal-degree factorisation: `pol` must be squarefree and a product of
irreducible polynomials of degree $d$, and $d$ must divide
$\deg(\mathit{pol})$. Appends the monic irreducible factors of `pol`,
each with exponent 1, to `factors`. This is `gr_poly_factor_equal_deg`;
see [](algorithms/polynomials/factoring-finite-fields#edf). An exception
is raised if $d$ does not divide the degree.
:::

## Irreducibility and squarefreeness {#irreducibility}

::: function
int fmpz_mod_poly_is_irreducible(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Returns whether `f` is irreducible over $\mathbb{F}_p$. As a special
case, constant polynomials, including zero, are considered irreducible.
This is `gr_poly_is_irreducible`, which uses the distinct-degree test for
degree below 600 and a filter for small factors followed by Rabin's test
for larger degrees; see [](algorithms/polynomials/irreducibility).
:::

::: function
int fmpz_mod_poly_is_irreducible_ddf(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_is_irreducible_rabin(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Return whether `f` is irreducible, using the distinct-degree test
(`gr_poly_is_irreducible_ddf`) or Rabin's test with the Frobenius powers
computed by modular composition (`gr_poly_is_irreducible_rabin`).
Constants are considered irreducible.
:::

::: function
int fmpz_mod_poly_is_irreducible_rabin_f(fmpz_t fac, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Rabin's irreducibility test for a possibly composite modulus $n$: either
sets `fac` to 1 and returns whether `f` is irreducible (when $n$ is
prime), or sets `fac` to a nontrivial factor of $n$, met as a
non-invertible leading coefficient, and returns an unspecified value.

::: details
The powers $x^{n^e} \bmod f$ are computed from a table of the powers
$x^{n^{2^i}}$ (`fmpz_mod_poly_frobenius_powers_2exp_precomp`). For `f`
of length at most 2 the function sets `fac` to 1 and returns 1. The
gcds are computed with `fmpz_mod_poly_gcd`, which can raise an exception
for a composite modulus.
:::
:::

::: function
int _fmpz_mod_poly_is_squarefree(const fmpz * f, slong len, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_is_squarefree(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Return whether `f` (given as `(f, len)` in the underscore version, with
any $\mathit{len} \ge 0$) is squarefree, by testing whether $\gcd(f, f')
= 1$ (`gr_poly_is_squarefree`). Nonzero constants are squarefree; as a
special case, the zero polynomial is not.
:::

::: function
int _fmpz_mod_poly_is_squarefree_f(fmpz_t fac, const fmpz * f, slong len, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_is_squarefree_f(fmpz_t fac, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Squarefreeness test for a possibly composite modulus $n$: either sets
`fac` to 1 and returns whether `f` is squarefree, or sets `fac` to a
nontrivial factor of $n$ and returns 0.
:::

## Roots {#roots}

::: function
void fmpz_mod_poly_roots(fmpz_mod_poly_factor_t r, const fmpz_mod_poly_t f, int with_multiplicity, const fmpz_mod_ctx_t ctx)

Sets `r` to the linear factors $x - a$ for the distinct roots $a \in
\mathbb{F}_p$ of the nonzero polynomial `f`. If `with_multiplicity` is
nonzero, the exponent of $x - a$ is the multiplicity of $a$ as a root of
`f`, otherwise it is 1. The modulus must be prime. A zero `f` causes an
exception. This is `gr_poly_roots_finite_field`; see
[](algorithms/polynomials/root-finding).
:::

::: function
int fmpz_mod_poly_roots_factored(fmpz_mod_poly_factor_t r, const fmpz_mod_poly_t f, int with_multiplicity, const fmpz_factor_t n, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_roots_factored_with_length_limit(fmpz_mod_poly_factor_t x0, const fmpz_mod_poly_t f, int with_mult, slong length_limit, const fmpz_factor_t fac, const fmpz_mod_ctx_t ctx)

Set the output to the linear factors $x - a$ for the distinct roots $a
\in \Z/n\Z$ of the nonzero polynomial `f`, where the `fmpz_factor_t`
holds the prime factorisation of the modulus $n$ (not checked). If the
multiplicity flag is nonzero, the exponent of $x - a$ is the largest $e$
such that $(x - a)^e$ divides $f$ in $(\Z/n\Z)[x]$, otherwise it is 1.
Return 1 on success and 0 if the roots were not computed because there
are too many of them. A zero `f` causes an exception.

::: details
The roots modulo each prime $p \mid n$ are found with
`fmpz_mod_poly_roots` (or are all residues if $f \equiv 0 \bmod p$),
lifted to all roots modulo the prime power $p^k \| n$ by Hensel lifting
with doubling precision, and combined by the Chinese remainder theorem;
see [](algorithms/polynomials/root-finding#prime-powers). The lists of
roots modulo a prime power are limited to about $2^{32}$ entries
($2^{25}$ on 32-bit systems), and the number of combined roots to the
same bound for `fmpz_mod_poly_roots_factored` and to `length_limit` for
`fmpz_mod_poly_roots_factored_with_length_limit`. If some prime power
gives no root, the functions return 1 with no roots, even if another
prime power has too many.
:::
:::
