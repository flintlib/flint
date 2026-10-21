---
title: Polynomial arithmetic over generic rings and its specialisations
kind: topic
summary: How the generic polynomial module gr_poly organises multiplication, division, Newton iteration, GCDs, evaluation, composition, power series and factorisation over any ring through overridable ring methods, and how nmod_poly, fmpz_poly, fmpz_mod_poly, fq_poly and arb_poly specialise or reuse these algorithms.
order: -10
---

Dense univariate polynomials appear in FLINT in two layers. The generic
module [`gr_poly`](gr_poly) implements the algorithms of univariate
polynomial arithmetic once, for coefficients in any ring given by a
`gr_ctx_t`: classical and fast multiplication, basecase and Newton
division, the Euclidean, half-gcd and subresultant algorithms,
multipoint evaluation and interpolation, composition, power series
functions and factorisation over finite fields. The type-specific
modules (`nmod_poly`, `fmpz_poly`, `fmpq_poly`, `fmpz_mod_poly`,
`fq_poly` and its variants, `arb_poly`, `acb_poly`, `ca_poly`) keep
hand-optimised kernels for the operations whose speed depends on the
coefficient representation, above all multiplication, and use the
generic implementations for much of the rest.

## Ring methods {#methods}

A handful of polynomial operations are *ring methods*: entries of the
method table of the coefficient ring
([](algorithms/generic/method-tables)), which a ring can override with
its own implementation. They are low products (`GR_METHOD_POLY_MULLOW`)
and middle products (`POLY_MULMID`), division with remainder
(`POLY_DIVREM`), power series inversion and division (`POLY_INV_SERIES`,
`POLY_DIV_SERIES`), square roots and exponentials of power series,
evaluation, Taylor shifts, GCDs and extended GCDs, resultants,
factorisation and roots, and the precomputation of moduli
(`POLY_PREINV_SET`). The default implementations of these methods are
the generic algorithms; the methods and the rings that override them
are listed in [](gr_poly#generic).

Every other function is written in terms of the ring methods and of
ring arithmetic. A ring that supplies a fast `POLY_MULLOW` therefore
gets fast Newton division, power series functions, subproduct trees,
composition and half-gcd from the generic code without further work,
and a ring that supplies its own GCD makes the squarefree
factorisation, the factorisation over finite fields and the root
finding of `gr_poly` use it. The element layout makes the two layers
interchangeable at the level of arrays: a `gr_poly` over the context
created by `gr_ctx_init_nmod` has exactly the coefficient array of an
`nmod_poly`, so wrappers between the layers are casts.

The generic code is written for rings in which not every operation
succeeds. Leading coefficients that are not provably nonzero are kept
(weak normalisation), divisions by non-invertible leading coefficients
fail with `GR_DOMAIN`, and undecidable zero tests over inexact rings
give `GR_UNABLE` ([](gr_poly#conventions),
[](algorithms/generic/inexact)). The type-specific modules assume an
exact ring with decidable equality and, over $\Z/n\Z$, usually a prime
modulus; they abort where the generic code returns a status.

## Multiplication {#multiplication}

The default multiplication method is the classical algorithm, written
so that each output coefficient is one dot product (`_gr_vec_dot_rev`).
With $O(n^2)$ ring operations this is the only algorithm that is valid
and efficient over every ring, and it automatically benefits from fast
dot products: delayed modular reduction for $\Z/n\Z$, exact
accumulation for $\Z$, a single rounding and a single error bound for
real balls. Squaring uses the symmetry of the product, and low and
middle products compute only the requested window
([](algorithms/polynomials/middle-product)).

`gr_poly` also provides one-level implementations of faster algorithms
that a ring can install in its method above a cutoff
([](gr_poly#mul-algorithms)): Karatsuba multiplication [@KO1962],
Toom-3 [@Too1963] with Bodrato's sequence of operations [@Bod2007], a
memory-efficient Toom-Cook algorithm of arbitrary order, Kronecker
substitution $y \mapsto x^k$ for bivariate products, which turns a
product in $R[y][x]$ into one product in $R[y]$ [@Kro1882], and the
reordering of a complex product into products of real and imaginary
parts. Their sub-products are done by the ring method, so they recurse
only through the ring's own dispatch. The coefficient rings override
the method as follows.

| Coefficient ring | Multiplication |
|---|---|
| `nmod` ($\Z/n\Z$, word-size $n$) | `nmod_poly` kernels: classical with modular dot products, Kronecker substitution (KS, KS2, KS4) and the small-prime FFT, chosen by length and modulus size ([](algorithms/polynomials/multiplication)) |
| `fmpz`, `fmpq` | `fmpz_poly` kernels: word-size loops, classical, Karatsuba, Kronecker substitution, Schönhage-Strassen and the small-prime FFT ([](algorithms/polynomials/multiplication-integers)); `fmpq_poly` clears denominators |
| `fmpz_mod` | product over $\Z$ with `fmpz_poly`, then reduction |
| `mpn_mod` (multi-word moduli) | classical or Karatsuba below a cutoff tabulated by modulus size, then the small-prime FFT, with Kronecker substitution as a fallback |
| `fq`, `fq_nmod`, `fq_zech` | classical below a cutoff, Kronecker substitution into an integer polynomial product above it (for `fq_zech`, conversion to the polynomial representation) |
| `arb`, `acb` | classical below 8 to 110 coefficients depending on the precision, otherwise the block algorithm, which splits the inputs into blocks of similar magnitude, multiplies them as exact integer polynomials, and bounds the error ([](algorithms/ball-arithmetic/polynomial-multiplication)) |
| `arf`, `acf` | classical below length 10, otherwise the ball algorithm applied to the floating-point inputs |
| `fmpzi` (Gaussian integers) | products of real and imaginary parts in $\Z[x]$, with Gauss's three-product trick for large coefficients |
| polynomial rings (`gr_poly` and `fmpz_poly` contexts) | bivariate Kronecker substitution when all lengths are at least 5 and the inner ring has a fast multiplication |
| `ca` | `ca_poly` multiplication |

Other rings, including number fields, `qqbar`, `nfloat` and the generic
quotient rings, use the classical algorithm. For rings with expensive
multiplication and cheap additions that lack an override, the
Karatsuba or Toom functions can be called directly.

::: summary algorithms/polynomials/multiplication :::

::: summary algorithms/polynomials/kronecker :::

## Division and Newton iteration {#division}

Division with remainder, power series inversion and division, square
roots, reciprocal square roots and exponentials of power series are
computed either by the classical (basecase) recurrence with $O(n^2)$
operations or by Newton iteration with $O(M(n))$ operations
[@Sie1972; @Kun1974; @BrentKung1978]. The generic code uses the
basecase for short inputs and Newton iteration above small fixed
cutoffs, since a Newton step is only as fast as the ring's
multiplication:

- `_gr_poly_divrem` uses the basecase when the divisor is a constant,
  when the quotient has length 1, or over finite rings when the
  quotient or the divisor has length at most 10, and otherwise computes
  the quotient by Newton iteration on the reversed polynomials and the
  remainder with one low product; `_gr_poly_div` uses the basecase when
  the quotient or the divisor has length at most 10, over any ring.
- `_gr_poly_inv_series` uses the basecase when the input has at most 8
  coefficients or the ring has only classical multiplication, and
  Newton iteration otherwise; `_gr_poly_div_series` uses the basecase
  when the divisor or the output has at most 8 coefficients (or the
  multiplication is classical) and otherwise Newton iteration with the
  multiplication by the dividend folded into the last step (the
  Karp-Markstein trick), using middle products where the ring provides
  them.
- Square roots and reciprocal square roots use a coupled Newton
  iteration, or J. C. P. Miller's recurrence when the input has at most
  8 coefficients and the output is longer, over an exact ring; the
  exponential uses the basecase recurrence when the input or the output
  has fewer than 10 terms (or the multiplication is classical) and Newton iteration with a basecase
  below length 30 otherwise.

Newton iteration needs an invertible leading coefficient of the divisor
(constant coefficient for power series). When it is not invertible, as
for a divisor $2x + 1$ over $\Z$, the generic division falls back to
basecase or divide-and-conquer algorithms that divide exactly by it at
each step, and succeed when the quotient has coefficients in the ring;
otherwise `GR_DOMAIN` is returned. Exact division (`gr_poly_divexact`)
assumes that the quotient exists and is the operation used by
fraction-free algorithms.

Precomputed moduli (`gr_poly_preinv_set`) store the inverse of the
reversed modulus or a sparse representation of it, and make the
repeated reductions in modular exponentiation, modular composition and
factorisation cheaper. The rings `nmod`, `fmpz_mod`, `mpn_mod`, the
`fq` family and `arb`/`acb` override the division methods: the finite
rings with their own basecase kernels and cutoffs around the generic
Newton iteration, and the balls with division through power series.

::: summary algorithms/polynomials/newton :::

## GCDs and resultants {#gcd}

The generic GCD chooses between the Euclidean algorithm, the half-gcd
algorithm and the subresultant pseudo-remainder sequence from the
properties of the ring: the half-gcd algorithm above length 200 over
finite fields, the Euclidean algorithm over other fields, and the
subresultant algorithm over unique factorisation domains such as $\Z$
and $\Z[y]$. Resultants follow the same pattern, with a division-free
Sylvester determinant as the last resort, and polynomial rings over
$\Z/p\Z$, $\Z$ and $\mathbb{Q}$ use bivariate multipoint and multimodular
algorithms. Over $\Z$, `fmpz_poly` has its own heuristic and
multimodular GCD ([](algorithms/polynomials/gcd-integers)) and
multimodular resultant ([](algorithms/polynomials/resultant-integers)).

::: summary algorithms/polynomials/gcd :::

## Evaluation and interpolation {#evaluation}

Single-point evaluation is a ring method whose default is Horner's
rule; `gr_poly` also offers rectangular splitting [@PS1973], which
replaces most full multiplications by scalar ones, a version that
evaluates with a few dot products (`gr_poly_evaluate_modular`), and
evaluation at an element of another ring, such as a matrix, or a ball
for an integer polynomial ([](gr_poly#evaluation)). Multipoint
evaluation and interpolation use the subproduct tree, with $O(M(n) \log
n)$ operations, or point-by-point Horner evaluation and Newton
interpolation for short inputs. `nmod_poly` and `fq_poly` call the
`gr_poly` subproduct-tree functions for fast multipoint evaluation.

::: summary algorithms/polynomials/subproduct-tree :::

## Composition and Taylor shifts {#composition}

Composition $f(g)$ uses Horner's rule for $f$ of length at most 7 and
the divide-and-conquer algorithm otherwise. The Taylor shift $f(x + c)$
is a ring method; the default uses Horner's scheme up to length 20 and
the divide-and-conquer algorithm above, which works over any ring.
`fmpz` and `arb`/`acb` override it with their own cutoffs, and the ball
version also uses the convolution algorithm, which needs division by
$n!$ and so characteristic zero.

Composition of power series $f(g) \bmod x^n$ uses Horner's rule for
short inputs, the divide-and-conquer algorithm when $f$ is short
relative to $n$ or the full composition is not much longer than $n$,
and otherwise the baby-step giant-step algorithm of Brent and Kung
[@BrentKung1978] below a cutoff (500 over exact rings, lower for balls
at low precision) and the algorithm of Kinoshita and Li [@KL2024] above
it. Series reversion uses Lagrange inversion with baby steps and giant
steps, or Newton iteration for long series. Modular composition $f(g)
\bmod h$ uses Horner's rule or the Brent-Kung algorithm, which reduces
to one matrix product over the ring and so profits from fast matrix
multiplication (`gr_mat_mul`). These algorithms are generic throughout:
`nmod_poly`, `fmpz_poly`, `fmpz_mod_poly`, `fq_poly` and `arb_poly`
call the `gr_poly` functions.

::: summary algorithms/polynomials/composition :::

## Power series functions {#series}

The elementary functions of power series (logarithm, exponential,
powers, trigonometric, hyperbolic and inverse trigonometric functions,
and Lambert's $W$) are built from the Newton iterations above and from
integration and differentiation, with basecase recurrences for short
series ([](gr_poly#special)). They work over any ring in which the
required divisions by integers and by the constant coefficient are
possible: $\mathbb{Q}$, finite fields of large characteristic (with
`GR_DOMAIN` when the characteristic divides a needed denominator), and
balls. Many `arb_poly` and `acb_poly` series functions are the
generic algorithms called with a ball context, and the special
functions of [`gr_series`](gr_series) are computed with them.

## Factorisation and roots {#factorisation}

Squarefree factorisation uses Yun's algorithm in characteristic zero
and the gcd iteration with $p$-th roots over finite fields. Over any
finite field whose context reports its order and characteristic (`nmod`
and `fmpz_mod` with prime moduli, `mpn_mod`, the `fq` family) `gr_poly`
provides a complete factorisation: distinct-degree factorisation by the
baby-step giant-step algorithm of Kaltofen and Shoup [@KalSho1998],
equal-degree factorisation by the Cantor-Zassenhaus algorithm
[@CanZas1981], irreducibility tests and root finding. These algorithms
rest on modular exponentiation and modular composition, hence on the
multiplication, division and matrix multiplication methods of the
ring. Factorisation and roots are ring methods: over $\Z$ and $\mathbb{Q}$
they use the Zassenhaus and van Hoeij algorithms of `fmpz_poly_factor`
([](algorithms/polynomials/factoring-integers), with
[](algorithms/polynomials/hensel-lifting)), and over the real and complex numbers the isolation of complex roots.
Conversely, `nmod_poly_factor` and `fmpz_mod_poly_factor` call the
generic finite-field algorithms.

::: summary algorithms/polynomials/squarefree :::

::: summary algorithms/polynomials/factoring-finite-fields :::

::: summary algorithms/polynomials/root-finding :::

## Specialised modules {#specialisations}

Each type-specific module combines its own kernels with the generic
algorithms.

| Module | Own kernels | Generic algorithms used |
|---|---|---|
| `nmod_poly` | multiplication and middle products, basecase division and power series kernels with delayed reduction, Euclidean GCD and resultant (also in Montgomery form), tuned cutoffs | Newton division and power series, half-gcd recursion, composition, power series composition and reversion, modular composition, fast multipoint evaluation, factorisation (through `nmod_poly_factor`) |
| `fmpz_poly`, `fmpq_poly` | multiplication, pseudo-division, heuristic and multimodular GCD, multimodular resultant, Hensel lifting and factorisation, real root isolation | composition, power series composition and reversion, power series inverse, powering |
| `fmpz_mod_poly` | multiplication through `fmpz_poly` | division, GCD and half-gcd, resultant, composition and modular composition, power series, factorisation |
| `fq_poly`, `fq_nmod_poly`, `fq_zech_poly` | multiplication by Kronecker substitution | division, GCD and half-gcd, composition and modular composition, evaluation, power series, factorisation |
| `arb_poly`, `acb_poly` | multiplication with error bounds, power series of special functions, root isolation | composition, Taylor shifts, reversion, many power series functions |

The functions of the specialised modules that call `gr_poly` are thin
wrappers: they create a context on the stack (for example with
`_gr_ctx_init_nmod`), pass the coefficient arrays, and call the generic
function, which in turn dispatches the core operations back to the
module's kernels through the ring methods.

When to use which layer:

- For a fixed coefficient type, the specialised module is usually at
  least as fast, needs no status checks, and has a richer interface (for
  example `fmpq_poly_gcd`, whose algorithm `gr_poly_gcd` over `fmpq`
  does not use).
- The generic layer is needed for coefficient rings without a module of
  their own (number fields, `qqbar`, `ca`, polynomial and series rings,
  matrices), for code that must work over several rings, and for rings
  where operations may fail (composite moduli, inexact numbers), where
  the status codes report the failure instead of aborting.

The guide [](gr_poly/getting-started) shows both kinds of use, and
[](gr_poly/gcd) the GCD over several rings.
