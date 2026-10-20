---
title: Multivariate polynomials: representations and algorithms
kind: topic
summary: An overview of FLINT's sparse multivariate polynomial modules: which coefficient rings they cover, how polynomials are represented, and which algorithms are used for multiplication, division, powering, square roots, GCDs, resultants and factorisation, with links to the detailed pages.
order: 0
---

FLINT represents multivariate polynomials sparsely, as sorted lists of
terms with packed exponent vectors, in one module per coefficient ring.
The algorithms are shared between the modules as far as possible: the
monomial code is in [](mpoly), and most GCD and factoring algorithms
exist in a version over $\Z$ and versions over finite fields, with the
rational and rational function modules built on top of the integer one.

| Coefficients | Polynomials | Factorisation | Fractions |
|---|---|---|---|
| any `gr` ring | `gr_mpoly` | | |
| $\Z$ | [](fmpz_mpoly) | [](fmpz_mpoly_factor) | [](fmpz_mpoly_q) |
| $\Q$ | [](fmpq_mpoly) (content times `fmpz_mpoly`) | [](fmpq_mpoly_factor) | |
| $\Z/n\Z$, word-size $n$ | `nmod_mpoly` | `nmod_mpoly_factor` | |
| $\Z/n\Z$, any $n$ | `fmpz_mod_mpoly` | `fmpz_mod_mpoly_factor` | `fmpz_mod_mpoly_q` |
| $\F_q$ | `fq_nmod_mpoly`, `fq_zech_mpoly` | `fq_nmod_mpoly_factor`, `fq_zech_mpoly_factor` | |

The generic `gr_mpoly` implements arithmetic (heap multiplication and
division, threaded versions) for any coefficient ring; the `gr` contexts
for the specialised types (`gr_ctx_init_fmpz_mpoly` and others) call the
specialised functions, which are faster and also provide GCDs and
factorisation.

## Representation {#representation}

::: summary algorithms/multivariate/representation :::

A polynomial is stored with its terms sorted for one of the orderings
`ORD_LEX`, `ORD_DEGLEX`, `ORD_DEGREVLEX`, and the exponent vectors packed
into words with a guard bit per field, so that comparing, multiplying and
testing divisibility of monomials are word operations. The packing width
grows automatically. The ordering affects the speed of some operations
(the dense array methods differ for lex and degree orderings) but not the
results of GCDs and factorisations, which are computed in lex order
internally and converted back.

## Arithmetic {#arithmetic}

::: summary algorithms/multivariate/multiplication :::

| Operation | Method | Used when | Functions |
|---|---|---|---|
| multiplication | Johnson's heap with chaining | sparse inputs, any ring | `fmpz_mpoly_mul_johnson`, `gr_mpoly_mul_heap` |
| | parallel heap | large inputs, several threads | `fmpz_mpoly_mul_heap_threaded` |
| | dense array accumulation | product dense in a box of moderate size | `fmpz_mpoly_mul_array` |
| | Kronecker substitution | dense inputs | `fmpz_mpoly_mul_dense` |
| powering | FPS recurrence | $\Z$, exponent $\ge 3$ | `fmpz_mpoly_pow_fps` |
| | repeated multiplication | finite fields | `nmod_mpoly_pow_rmul` |

`fmpz_mpoly_mul` and `nmod_mpoly_mul` choose automatically; the rules are
described in [](fmpz_mpoly/mul).

::: summary algorithms/multivariate/division :::

| Operation | Method | Functions |
|---|---|---|
| exact division | heap division with early termination; threads for large dividends | `fmpz_mpoly_divides` |
| division with remainder | heap division (floor division of coefficients over $\Z$) | `fmpz_mpoly_divrem` |
| pseudo-division | heap division with scaling | `fmpz_mpoly_quasidivrem` |
| division by several polynomials | heap | `fmpz_mpoly_divrem_ideal` |
| square root | heap, with non-square certificates | `fmpz_mpoly_sqrt` |

## GCDs and resultants {#gcd}

::: summary algorithms/multivariate/gcd :::

`fmpz_mpoly_gcd` first reduces the problem (monomial contents, deflation,
missing variables, univariate cases, divisibility tests) and then chooses
between Brown's dense modular algorithm (two variables, dense inputs),
the Berlekamp–Massey based sparse algorithm of Hu and Monagan (sparse
inputs in three or more variables) and Zippel's sparse algorithm, using
degree and term count estimates from random images. Hensel lifting
(EEZ-GCD) and subresultants are available on request. Over finite
fields, `nmod_mpoly_gcd` also uses Hensel lifting, and computes over
extension fields when the field is too small:

::: summary algorithms/multivariate/modp-gcd :::

::: summary algorithms/multivariate/sparse-interpolation :::

Resultants and discriminants (`fmpz_mpoly_resultant`,
`fmpz_mpoly_discriminant`) use the subresultant algorithm of Ducos
[@Duc2000] on univariate polynomials with multivariate coefficients, except
that bivariate dense inputs go through a multimodular evaluation and
interpolation method:

::: summary algorithms/multivariate/modp-bivariate-resultant :::

## Factorisation {#factorisation}

::: summary algorithms/multivariate/factorisation :::

::: summary algorithms/multivariate/hensel-lifting :::

| Stage | Method | Functions |
|---|---|---|
| content | GCDs of coefficients with respect to each variable | `fmpz_mpoly_factor_content` |
| squarefree | Yun's algorithm (with $p$-th roots in characteristic $p$) | `fmpz_mpoly_factor_squarefree` |
| simplification | unimodular compression of exponents; Newton polytope test | `mpoly_compression_set` |
| low degree | linear and quadratic in the main variable | `fmpz_mpoly_factor` |
| bivariate dense | Hensel lifting and recombination | `fmpz_bpoly_factor` |
| general | evaluation, degree pruning, then Wang (dense lifting), Zippel (sparse lifting) or Zassenhaus (lifting with recombination) | `fmpz_mpoly_factor_wang`, `fmpz_mpoly_factor_zippel`, `fmpz_mpoly_factor_zassenhaus` |

Factoring over finite fields follows the same stages, with the changes
described in:

::: summary algorithms/multivariate/modp-factoring :::

## Practical advice {#advice}

- Use the specialised modules rather than `gr_mpoly` over $\Z$, $\Q$ and
  finite fields: they accumulate coefficients without
  intermediate reductions, use dense methods when they pay, and provide
  GCDs and factorisation.
- The ordering only matters for the leading terms (and hence for division
  with remainder and Gröbner bases) and for the speed of the array
  methods; lex order is the natural one for recursive computations.
- Threads (`flint_set_num_threads`) are used by multiplication and exact
  division of large polynomials, and by Brown's GCD algorithm in
  `nmod_mpoly_gcd`.
- GCDs and factorisations are randomised but never return wrong results:
  every algorithm verifies its output, by exact division or by a proof
  of its stopping criterion. Their running time can vary with the random
  choices.
