---
title: Linear algebra over generic rings and its specialisations
kind: topic
summary: How the generic matrix module gr_mat organises multiplication, LU and fraction-free LU decomposition, solving, determinants, ranks, echelon forms, characteristic polynomials and eigenvalues over any ring, how it handles rings that are not fields and inexact rings, and how nmod_mat, fmpz_mat, mpn_mod, arb_mat, acb_mat and fq_mat specialise or reuse these algorithms.
order: -10
---

Dense linear algebra in FLINT comes in two layers, as polynomial
arithmetic does ([](algorithms/polynomials/overview)). The generic
module [`gr_mat`](gr_mat) implements matrix multiplication, Gaussian
elimination and its fraction-free variant, triangular solving,
determinants, ranks, echelon forms, nullspaces, inverses,
characteristic and minimal polynomials and eigenvalues once, for
entries in any ring given by a `gr_ctx_t`. The type-specific modules
(`nmod_mat`, `fmpz_mat`, `fmpq_mat`, `fmpz_mod_mat`, `fq_mat` and its
variants, `arb_mat`, `acb_mat`, `ca_mat`) provide kernels whose speed
depends on the entry representation, above all matrix multiplication,
and algorithms that exploit the structure of their ring, such as
multimodular and $p$-adic methods over $\Z$ and preconditioned
elimination for balls.

## Ring methods and the role of multiplication {#methods}

A few matrix operations are ring methods that a coefficient ring can
override ([](algorithms/generic/method-tables)): multiplication
(`GR_METHOD_MAT_MUL`), the LU decomposition, the determinant, the
characteristic polynomial, pivot search, row reduction against an
echelon basis (`MAT_REDUCE_ROW`, used by the minimal polynomial),
triangular solving, the exponential and diagonalisation. The rest of
`gr_mat` is written in terms of these and of ring arithmetic.

Matrix multiplication matters most. The recursive LU decomposition,
triangular solving with many right-hand sides, inversion and the
algorithms built on them do their work in matrix products, so a ring
that overrides `MAT_MUL` with a fast algorithm makes them run in
$O(n^\omega)$ operations with the exponent of its multiplication. The
default is the classical algorithm, with each entry one dot product
(`_gr_vec_dot`); it thus benefits from fast dot products (delayed
modular reduction, exact accumulation, one rounding per entry for
floating-point and ball arithmetic). One-level Strassen-Winograd,
Waksman and Rosowski multiplications are available for rings to install
above a cutoff [@Str1969; @Win1971; @Waksman1970; @Rosowski2023]. The
overrides by ring, from `nmod_mat_mul` over $\Z/n\Z$ to the block
algorithm of `arb_mat_mul`, are tabulated in [](gr_mat/mul).

::: summary algorithms/linear-algebra/matrix-multiplication :::

::: summary algorithms/linear-algebra/strassen-winograd :::

## Gaussian elimination and LU decomposition {#lu}

Over a field, `gr_mat_lu` computes a decomposition $PA = LU$ with a row
permutation $P$, unit lower triangular $L$ and upper triangular $U$
(in reduced row echelon shape for singular or rectangular matrices).
The default `gr_mat_lu_generic` uses classical elimination when the
matrix has fewer than 5 rows or columns and otherwise the recursive
algorithm, which splits the columns in two halves, decomposes the left
half, updates the right half with one triangular solve and one matrix
product, and recurses; with fast multiplication its cost is
$O(n^\omega)$ ([](algorithms/linear-algebra/gaussian-elimination)).

Pivots are chosen by the ring method `MAT_FIND_NONZERO_PIVOT`. The
default takes the first entry that is provably nonzero, and if there
is none returns `GR_UNABLE` when some entry cannot be decided and
`GR_DOMAIN` otherwise; the floating-point
rings (`arf`, `acf`, `nfloat`) choose the entry of largest magnitude
for numerical stability, and `ca` works harder to decide zero tests.
The rings `fmpz_mod`, `mpn_mod`, `nfloat` and `nfloat_complex`
override the LU decomposition itself with their own cutoffs and
kernels.

Over an integral domain that is not a field, the fraction-free LU
decomposition (`gr_mat_fflu`) divides only exactly, using Sylvester's
identity as in Bareiss' algorithm [@Bar1968]: all intermediate entries
are minors of the input, so their size is bounded and no fractions
appear ([](algorithms/linear-algebra/fraction-free)).

The generic functions built on these decompositions choose between
them by the properties of the ring (`gr_mat_nonsingular_solve` and
`gr_mat_inv` are ring methods, which a ring can replace; the $p$-adic
`padic_radix` rings use elimination with full pivoting by valuation
levels for the optimal precision):

| Operation | Algorithm |
|---|---|
| `gr_mat_nonsingular_solve` | LU over fields (and rings pretending to be fields) and over inexact rings, fraction-free LU otherwise (which succeeds only if the solution has entries in the ring) |
| `gr_mat_nonsingular_solve_den` | fraction-free LU: $X$ and $d$ with $AX = dB$ |
| `gr_mat_solve_field` | LU with a reduction of the right-hand side; general (rectangular or singular) systems over a field |
| `gr_mat_inv` | `gr_mat_nonsingular_solve` with the identity matrix as right-hand side |
| `gr_mat_rank` | LU over finite fields and finite rings pretending to be fields, fraction-free LU over other integral domains, `GR_DOMAIN` over rings with zero divisors |
| `gr_mat_rref`, `gr_mat_rref_den` | LU, fraction-free LU |
| `gr_mat_nullspace` | from the reduced echelon form, fraction-free over integral domains that are not fields |

::: summary algorithms/linear-algebra/gaussian-elimination :::

::: summary algorithms/linear-algebra/fraction-free :::

## Determinants {#det}

The determinant is a ring method. The generic choices are cofactor
expansion up to $4 \times 4$, and then:

- over any commutative ring, Berkowitz's division-free algorithm
  [@Ber1984], with $O(n^4)$ operations (`gr_mat_det_generic`, the
  default);
- over fields, the LU decomposition, with $O(n^\omega)$ operations
  (`gr_mat_det_generic_field`, used by the floating-point rings);
- over integral domains, the fraction-free LU decomposition, with
  $O(n^3)$ exact divisions (`gr_mat_det_generic_integral_domain`).

The integers and rationals use the multimodular algorithm of
`fmpz_mat_det` ([](algorithms/linear-algebra/multimodular)), `nmod`
calls `nmod_mat_det`, `nmod8`, `nmod32`, `fmpz_mod` and `mpn_mod` use LU
with Berkowitz's algorithm as the fallback for composite moduli, the
`fq` rings use LU, and balls use Gaussian elimination in ball arithmetic
with approximate preconditioning for large matrices. Because the
default is division-free, the determinant over a field without an
override costs $O(n^4)$ operations; [](gr_mat/det) gives the details.
The permanent, by
contrast, has no polynomial-time algorithm [@Val1979]; `gr_mat_permanent`
uses Ryser's and Glynn's formulas [@Rys1963; @Gly2010].

## Characteristic and minimal polynomials {#charpoly}

The default characteristic polynomial is Berkowitz's algorithm, which
is division-free and works over any commutative ring. `gr_mat` also
implements Danilevsky's reduction to companion form [@Dan1937], with
$O(n^3)$ operations over a field, the Faddeev-LeVerrier algorithm
[@FadSom1972] and its baby-step giant-step variant, which need division
by $1, \ldots, n$, and reduction to Hessenberg form by Gaussian or
Householder transformations followed by the Hessenberg recurrence
[@Wil1965]. The rings override the method with their own choice:
`nmod` uses Danilevsky's algorithm above dimension 8 and `fmpz_mod`
from dimension 16, both with Berkowitz's algorithm as the fallback when
a pivot is not invertible (composite moduli); the `fq` rings use
Danilevsky's algorithm above dimension 12; and `fmpz` and `fmpq` use
`fmpz_mat_charpoly`, which is Berkowitz's algorithm for small matrices
and a multimodular algorithm otherwise.

The minimal polynomial over a field (`gr_mat_minpoly_field`) is computed
from Krylov sequences $v, Av, A^2v, \ldots$ for the unit vectors $v$,
reducing each new vector against an echelon basis of the previous ones
(the ring method `MAT_REDUCE_ROW`) to find the first linear relation,
and taking the least common multiple of the relations.

Eigenvalues are the roots of the characteristic polynomial, computed
with `gr_poly_roots` in the ring itself or in another ring
(`gr_mat_eigenvalues_other`, for example the complex algebraic numbers
for a rational matrix). For complex balls, `acb` overrides the
diagonalisation method: an approximate eigendecomposition by the QR
algorithm is certified with `acb_mat_eig_simple`, which succeeds when
the eigenvalues are simple and well separated.

::: summary algorithms/linear-algebra/charpoly :::

## Rings that are not fields, and inexact rings {#rings}

The generic functions never return a wrong answer when the ring lacks
a property they need; they return a status instead. Over a ring with
zero divisors such as $\Z/6\Z$, the elimination-based functions fail
with `GR_DOMAIN` (or `GR_UNABLE`) when they meet a non-invertible pivot,
while the division-free determinant, characteristic polynomial and
adjugate still succeed. A ring $\Z/n\Z$ that pretends to be a field runs
the field algorithms and reports a factor of $n$ if it finds one.

Over real and complex balls, a pivot whose ball contains zero cannot be
used and cannot be excluded, so rank, solving and inversion return
`GR_UNABLE` when the outcome depends on it, while an exactly zero
pivot is recognised. Division-free algorithms and ball determinants
always return an enclosure. Over floating-point rings (`arf`, `acf`,
`nfloat`), values are treated as exact and pivoting by magnitude gives
the usual numerical algorithms, without error bounds. The guide
[](gr_mat/exactness) works through these cases with examples.

## Specialised modules {#specialisations}

| Module | Own algorithms | Generic algorithms used |
|---|---|---|
| `nmod_mat` | multiplication (classical with delayed reduction, SIMD and floating-point kernels, Strassen-Winograd), LU decomposition (classical, delayed, recursive), solving, inverse, rank, echelon form, nullspace, determinant (Howell form for composite moduli), Danilevsky and Berkowitz characteristic polynomial, minimal polynomial | characteristic polynomial dispatch |
| `fmpz_mat`, `fmpq_mat` | multiplication (word-size kernels, multimodular, Waksman, Strassen, small-prime FFT), multimodular determinant, rank and characteristic polynomial, $p$-adic (Dixon) solving, fraction-free echelon forms, Hermite and Smith normal forms, LLL | Berkowitz characteristic polynomial for small inputs, permanent, simple operations |
| `fmpz_mod_mat` | multiplication through `fmpz_mat` | LU, solving, rank, echelon form, nullspace, determinant, characteristic and minimal polynomials |
| `fq_mat`, `fq_nmod_mat`, `fq_zech_mat` | multiplication by Kronecker substitution | LU, solving, inverse, rank, echelon form, nullspace, characteristic and minimal polynomials |
| `mpn_mod` | multiplication (classical, Waksman, multimodular), LU, determinant, row reduction | the rest of `gr_mat` |
| `arb_mat`, `acb_mat` | block multiplication with exact integer products, LU and preconditioned solving and determinants, inverse, Cholesky and $LDL^T$ decompositions, exponential, certified eigenvalues | characteristic polynomial (Berkowitz), companion matrices, powers, special matrices |
| `ca_mat` | exact solving with zero tests over number fields, Jordan forms | most algorithms, called with the `ca` context |

The specialised algorithms are described in
[](algorithms/linear-algebra/multimodular) (integers and rationals),
[](algorithms/linear-algebra/integer-matrix-multiplication),
[](algorithms/linear-algebra/small-modulus-matrix-multiplication),
[](algorithms/linear-algebra/lll),
[](algorithms/linear-algebra/hermite-smith),
[](algorithms/linear-algebra/polynomial-matrices) and
[](algorithms/exact/linear-algebra) (exact real and complex numbers).
Over $\Z$ and $\mathbb{Q}$, the ring methods make `gr_mat_mul`,
`gr_mat_det`, `gr_mat_charpoly` and `gr_mat_permanent` call
`fmpz_mat_mul`, `fmpz_mat_det`, `_fmpz_mat_charpoly` and
`fmpz_mat_permanent` (or their `fmpq_mat` counterparts).

Some specialised algorithms are not reached through `gr_mat`: over
`arb` and `acb`, `gr_mat_nonsingular_solve` and `gr_mat_inv` use the
generic ball LU decomposition (as `arb_mat_solve_lu` does) rather than
the preconditioned solver of `arb_mat_solve`, which gives much tighter
enclosures for large well-conditioned systems; and over `fmpz`,
`gr_mat_nonsingular_solve_den` uses fraction-free elimination rather
than the $p$-adic solver of `fmpz_mat_solve`. For a fixed entry type and large
inputs, the specialised module is then the better choice; the generic
layer is needed for entry rings without a module of their own
(polynomial rings, number fields, `qqbar`, `ca`), for code that works
over several rings, and where the status codes are needed.
