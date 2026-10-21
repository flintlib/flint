---
title: nmod_poly_factor.h
header: nmod_poly_factor.h
description: Factorisation, irreducibility and squarefreeness testing, and root finding for polynomials over $\Z/n\Z$ with word-size modulus.
category: Univariate polynomials
order: 6
---

This module factors polynomials over the prime field $\mathbb{F}_p$, $p <
2^{\bits}$, represented as `nmod_poly_t`, tests them for irreducibility
and squarefreeness, and finds their roots, over $\mathbb{F}_p$ and also
over $\Z/n\Z$ for composite $n$ with known factorisation. A factorisation
is stored in an `nmod_poly_factor_t`, a list of polynomials with
exponents.

The factorisation, irreducibility and root-finding functions are thin
wrappers around the generic finite-field implementations in `gr_poly`
(`gr_poly_factor_finite_field`, `gr_poly_is_irreducible`,
`gr_poly_roots_finite_field` and the functions they call), which view
the `nmod_poly_t` as a `gr_poly_t` over the `nmod` context without
copying it. A complete factorisation consists of a squarefree
factorisation followed, for each squarefree part, by either the
Cantor–Zassenhaus algorithm with classical distinct-degree factorisation
(small degrees and small $p$) or the Kaltofen–Shoup algorithm with
baby-step giant-step distinct-degree factorisation and modular
composition; Berlekamp's algorithm is available on request. The guide
[](nmod_poly_factor/factor) explains the choices and their costs with
examples, and the algorithms are described in detail in:

::: summary algorithms/polynomials/factoring-finite-fields :::

::: summary algorithms/polynomials/root-finding :::

## Types and conventions {#conventions}

::: type
nmod_poly_factor_struct
nmod_poly_factor_t

A polynomial in factored form. The structure has the fields `p` (an
array of `alloc` polynomials, of which the first `num` are used), `exp`
(an array of `alloc` exponents of type `slong`), `num` and `alloc`.
Entry $i$ for $0 \le i < \mathit{num}$ stands for the factor
$\mathtt{p[i]}^{\mathtt{exp[i]}}$. An `nmod_poly_factor_t` is an array of
one `nmod_poly_factor_struct`, so it is passed by reference.
:::

- The factorisation, irreducibility, squarefreeness and root-finding
  functions require the modulus $p$ to be prime, except
  `nmod_poly_roots_factored`. This is not checked: a composite modulus
  gives wrong results or an exception.
- The factors computed by this module are monic. `nmod_poly_factor` and
  its variants return the leading coefficient of the input separately.
- The factors are stored in no particular order. The algorithms are
  randomised, but every function seeds its own random state with a fixed
  value, so the output is reproducible.
- Most functions overwrite the output factorisation, but
  `nmod_poly_factor_cantor_zassenhaus`, `nmod_poly_factor_berlekamp`,
  `nmod_poly_factor_kaltofen_shoup`, `nmod_poly_factor_distinct_deg` and
  `nmod_poly_factor_equal_deg` append to it; each function says which.
- Distinct-degree and equal-degree factorisation of polynomials of
  degree above about 200 use several threads when they have been made
  available with `flint_set_num_threads`.

## Memory management {#memory}

::: function
void nmod_poly_factor_init(nmod_poly_factor_t fac)

Initialises `fac` as the empty factorisation ($\mathit{num} = 0$), with
space for 5 factors.
:::

::: function
void nmod_poly_factor_clear(nmod_poly_factor_t fac)

Frees the memory used by `fac`, including all allocated polynomials.
:::

::: function
void nmod_poly_factor_realloc(nmod_poly_factor_t fac, slong alloc)

Changes the space of `fac` to exactly `alloc` factors. For `alloc` = 0
the factorisation is cleared and reinitialised (with space for 5 factors
and $\mathit{num} = 0$). If `alloc` is smaller than the number of
factors, the factors beyond the first `alloc` are discarded.
:::

::: function
void nmod_poly_factor_fit_length(nmod_poly_factor_t fac, slong len)

Ensures that `fac` has space for at least `len` factors. If it has to
grow, it grows to at least twice its previous allocation.
:::

::: function
void nmod_poly_factor_swap(nmod_poly_factor_t a, nmod_poly_factor_t b)

Swaps `a` and `b` efficiently.
:::

## Basic manipulation {#manipulation}

::: function
void nmod_poly_factor_set(nmod_poly_factor_t res, const nmod_poly_factor_t fac)

Sets `res` to a copy of `fac`, including the moduli of the factors.
:::

::: function
void nmod_poly_factor_get_poly(nmod_poly_t a, const nmod_poly_factor_t b, slong i)

Sets `a` to the factor `b->p[i]`, for $0 \le i < \mathit{num}$. The
exponent is `b->exp[i]`.
:::

::: function
void nmod_poly_factor_insert(nmod_poly_factor_t fac, const nmod_poly_t poly, slong exp)

Multiplies the factorisation by $\mathit{poly}^{\mathit{exp}}$: if
`poly` is equal to one of the factors, `exp` is added to its exponent,
otherwise a copy of `poly` is appended with exponent `exp`. Constant
polynomials are ignored. The search is linear in the number of factors.
:::

::: function
void nmod_poly_factor_concat(nmod_poly_factor_t res, const nmod_poly_factor_t fac)

Multiplies `res` by the factorisation `fac`, by inserting each factor of
`fac` with `nmod_poly_factor_insert`. `res` and `fac` must be distinct
objects.
:::

::: function
void nmod_poly_factor_pow(nmod_poly_factor_t fac, slong exp)

Multiplies every exponent of `fac` by `exp`, which raises the factored
polynomial to the power `exp`.
:::

## Input and output {#io}

::: function
void nmod_poly_factor_print(const nmod_poly_factor_t fac)
void nmod_poly_factor_print_pretty(const nmod_poly_factor_t fac, const char * var)

Print the factors to standard output, one per line, as the polynomial
followed by ` ^ ` and the exponent. `nmod_poly_factor_print` writes the
polynomials in the format of `nmod_poly_print`, and
`nmod_poly_factor_print_pretty` in the format of `nmod_poly_print_pretty`
with the variable name `var`.
:::

## Factorisation {#factorisation}

::: function
ulong nmod_poly_factor(nmod_poly_factor_t result, const nmod_poly_t input)

Sets `result` to the factorisation of `input` into monic irreducible
polynomials over $\mathbb{F}_p$ and returns the leading coefficient of
`input`, so that `input` equals the leading coefficient times the
product of the factors raised to their exponents. A constant `input`
gives the empty factorisation and returns the constant (0 for the zero
polynomial).

::: details
This is `gr_poly_factor_finite_field` over the `nmod` context. If `input`
is a polynomial in $x^m$ for some $m \ge 2$, it is first deflated, and
the factors $g(x)$ of the deflated polynomial are replaced by the factors
of $g(x^m)$. The squarefree factors of the (deflated) input of degree
$n$ are then split by the Cantor–Zassenhaus algorithm if $n <
125/(\log_2 p)^{1.3}$ (for example $n < 125$ for $p = 2$ and $n \le 32$
for $p = 7$), and by the Kaltofen–Shoup algorithm otherwise, which means
always for primes of 40 bits or more. See [](nmod_poly_factor/factor).
:::
:::

::: function
ulong nmod_poly_factor_with_cantor_zassenhaus(nmod_poly_factor_t result, const nmod_poly_t input)
ulong nmod_poly_factor_with_berlekamp(nmod_poly_factor_t result, const nmod_poly_t input)
ulong nmod_poly_factor_with_kaltofen_shoup(nmod_poly_factor_t result, const nmod_poly_t input)

Like `nmod_poly_factor`, but split the squarefree factors with the
Cantor–Zassenhaus algorithm, Berlekamp's algorithm or the Kaltofen–Shoup
algorithm respectively. Deflation is still applied.
:::

::: function
void nmod_poly_factor_cantor_zassenhaus(nmod_poly_factor_t res, const nmod_poly_t f)
void nmod_poly_factor_berlekamp(nmod_poly_factor_t res, const nmod_poly_t f)
void nmod_poly_factor_kaltofen_shoup(nmod_poly_factor_t res, const nmod_poly_t f)

Append the factorisation of `f` into monic irreducible polynomials, with
exponents, to `res`, using the Cantor–Zassenhaus algorithm, Berlekamp's
algorithm or the Kaltofen–Shoup algorithm (with deflation, as in
`nmod_poly_factor`). The leading coefficient of `f` is discarded;
nothing is appended for a constant `f`.
:::

::: function
void nmod_poly_factor_squarefree(nmod_poly_factor_t res, const nmod_poly_t f)

Sets `res` to the squarefree factorisation of `f`: monic squarefree
polynomials $g_i$, pairwise coprime and with distinct exponents $e_i$,
such that $f = c \prod_i g_i^{e_i}$ with $c$ the leading coefficient of
$f$. A constant `f` gives the empty factorisation.

::: details
This is `gr_poly_factor_squarefree`: Yun's algorithm for the part with
multiplicities prime to $p$, and recursion on the $p$-th root of the
remaining part, which is a polynomial in $x^p$. See
[](algorithms/polynomials/squarefree).
:::
:::

::: function
void nmod_poly_factor_distinct_deg(nmod_poly_factor_t res, const nmod_poly_t poly, slong * const * degs)
void nmod_poly_factor_distinct_deg_threaded(nmod_poly_factor_t res, const nmod_poly_t poly, slong * const * degs)

Distinct-degree factorisation of a squarefree polynomial `poly` of
degree at least 1, which need not be monic. Appends to `res` monic
polynomials $f_1, \ldots, f_k$, each with exponent 1, such that $f_j$ is
the product of all monic irreducible factors of `poly` of one degree
$d_j$, and writes $d_j$ to `(*degs)[j - 1]`. The array `*degs` must have
room for $\deg(\mathit{poly})$ entries; the $d_j$ are distinct but not
sorted. If `poly` is not squarefree the output is meaningless.

::: details
This is `gr_poly_factor_distinct_deg`, the baby-step giant-step
algorithm of Kaltofen and Shoup with Brent–Kung modular composition
[@KalSho1998]; see
[](algorithms/polynomials/factoring-finite-fields#ddf). It uses several
threads by itself when they are available, so
`nmod_poly_factor_distinct_deg_threaded` is now an alias of
`nmod_poly_factor_distinct_deg`. A constant `poly` causes an exception.
:::
:::

::: function
void nmod_poly_factor_equal_deg(nmod_poly_factor_t factors, const nmod_poly_t pol, slong d)

Equal-degree factorisation: `pol` must be squarefree and a product of
irreducible polynomials of degree $d$, and $d$ must divide
$\deg(\mathit{pol})$. Appends the monic irreducible factors of `pol`,
each with exponent 1, to `factors`.

::: details
This is `gr_poly_factor_equal_deg`: Cantor–Zassenhaus splitting of
random elements of the Berlekamp subalgebra obtained with the trace
map, which is computed by modular composition, and for many factors of
large degree the splitting by the minimal polynomial of one such element;
see [](algorithms/polynomials/factoring-finite-fields#edf). An exception
is raised if $d$ does not divide the degree. If `pol` does not have the
required form, the function can fail with an exception or return
reducible factors.
:::
:::

## Irreducibility and squarefreeness {#irreducibility}

::: function
int nmod_poly_is_irreducible(const nmod_poly_t f)

Returns whether `f` is irreducible over $\mathbb{F}_p$. As a special
case, constant polynomials, including zero, are considered irreducible
(unlike `gr_poly_is_irreducible`, which returns false for them).

::: details
For $p \le \max(200, 2\,\mathrm{len}(f))$ the function first looks for a
root by evaluating $f$ at every element of $\mathbb{F}_p$. For degree
below 600 it then uses `nmod_poly_is_irreducible_ddf`; for higher
degrees it checks squarefreeness, looks for factors of degree at most 4
with $\gcd(f, x^{p^i} - x)$, and finishes with Rabin's test computed by
modular composition (`gr_poly_is_irreducible_rabin`). See
[](algorithms/polynomials/irreducibility).
:::
:::

::: function
int nmod_poly_is_irreducible_ddf(const nmod_poly_t f)

Returns whether `f` is irreducible, by the coarse phase of the
baby-step giant-step distinct-degree factorisation: after a
squarefreeness test, $f$ is irreducible if and only if it has no factor
of degree at most $\deg(f)/2$. Constants are considered irreducible.
This is `gr_poly_is_irreducible_ddf`; it stops at the first interval of
degrees that contains a factor.
:::

::: function
int nmod_poly_is_irreducible_rabin(const nmod_poly_t f)

Returns whether `f` is irreducible, by Rabin's test [@Rab1980b]: $f$ of
degree $n$ is irreducible if and only if $x^{p^n} \equiv x \pmod f$ and
$\gcd(x^{p^{n/r}} - x, f) = 1$ for every prime $r \mid n$. Constants are
considered irreducible.

::: details
Unlike the other irreducibility tests, this function is not a wrapper of
the generic implementation: it computes each power $x^{p^e} \bmod f$ by
$e$ successive exponentiations by $p$, at a cost of about $e \log_2 p$
multiplications modulo $f$, which makes it much slower than
`nmod_poly_is_irreducible` for large degrees.
:::
:::

::: function
int _nmod_poly_is_squarefree(nn_srcptr f, slong len, nmod_t mod)
int nmod_poly_is_squarefree(const nmod_poly_t f)

Return whether `f` (given as `(f, len)` in the underscore version, with
any $\mathit{len} \ge 0$) is squarefree, by testing whether $\gcd(f, f')
= 1$; if $f' = 0$ and $\deg f \ge 2$, then $f$ is a $p$-th power and not
squarefree. Nonzero constants are squarefree; as a special case, the
zero polynomial is not.
:::

## Roots {#roots}

::: function
void nmod_poly_roots(nmod_poly_factor_t r, const nmod_poly_t f, int with_multiplicity)

Sets `r` to the linear factors $x - a$ for the distinct roots $a \in
\mathbb{F}_p$ of the nonzero polynomial `f`. If `with_multiplicity` is
nonzero, the exponent of $x - a$ is the multiplicity of $a$ as a root of
`f`, otherwise it is 1. The modulus must be prime. A zero `f` causes an
exception.

::: details
This is `gr_poly_roots_finite_field`: squarefree factorisation, then for
each squarefree part $g$ the split $\gcd(g, x^{(p-1)/2} \mp 1)$ into the
roots that are squares and those that are non-squares, followed by
equal-degree factorisation with $d = 1$. See
[](algorithms/polynomials/root-finding) and [](nmod_poly_factor/roots).
:::
:::

::: function
int nmod_poly_roots_factored(nmod_poly_factor_t r, const nmod_poly_t f, int with_multiplicity, const n_factor_t * n)

Sets `r` to the linear factors $x - a$ for the distinct roots $a \in
\Z/n\Z$ of the nonzero polynomial `f`, where `n` holds the prime
factorisation of the modulus $n$ of `f` (not checked). If
`with_multiplicity` is nonzero, the exponent of $x - a$ is the largest
$e$ such that $(x - a)^e$ divides $f$ in $(\Z/n\Z)[x]$, otherwise it is
1. Returns 1 on success and 0 if the roots were not computed because
there are too many of them (more than about $2^{32}$, or $2^{25}$ on
32-bit systems). A zero `f` causes an exception.

::: details
The roots modulo each prime $p \mid n$ are found with `nmod_poly_roots`
(or are all of $\mathbb{F}_p$ if $f \equiv 0 \bmod p$), lifted to roots
modulo the full prime power $p^k \| n$ by Hensel lifting, where each root
modulo $p^e$ can lift to several roots modulo $p^{2e}$ or to none, and
combined by the Chinese remainder theorem. If some prime power gives no
root, the function returns 1 with no roots, even if another prime power
has too many. See [](algorithms/polynomials/root-finding#prime-powers).
:::
:::
