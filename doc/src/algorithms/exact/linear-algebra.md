---
title: Linear algebra over exact numbers
summary: Gaussian elimination, fraction-free elimination, division-free determinants and characteristic polynomials, eigenvalues, Jordan forms and matrix functions over exact real and complex numbers, where every pivot decision is a zero test that may fail; how ca_mat builds on the generic gr_mat algorithms and specialises to rational and number field matrices.
complexity: Elimination costs $O(n^3)$ field operations and up to $O(n^2)$ zero tests; the Berkowitz algorithm costs $O(n^4)$ ring operations and no zero tests; eigenvalue and Jordan form computations add root finding for the characteristic polynomial and $O(n)$ rank computations per eigenvalue.
order: 60
---

Linear algebra over a field needs to decide, at each step of Gaussian
elimination, whether a candidate pivot is zero. Over exact real and
complex numbers ([](ca)) this decision is a zero test that may fail
(see [](algorithms/exact/zero-testing)), and arithmetic is expensive and
grows the representations. The [](ca_mat) module therefore relies on
the generic algorithms of [](gr_mat), which report failure through their
status codes, combined with a careful choice of pivots, algorithms that
avoid zero tests, and conversions to faster types when the entries
permit.

## Generic algorithms {#generic}

Most functions of [](ca_mat) are wrappers that call the corresponding
`gr_mat` function with a `ca` context: `ca_mat_lu` calls
`gr_mat_lu_generic` (classical Gaussian elimination below dimension 5,
block recursive elimination above), `ca_mat_rank`, `ca_mat_rref`,
`ca_mat_right_kernel`, `ca_mat_nonsingular_solve_lu`, the Berkowitz and
Danilevsky characteristic polynomials, the Jordan form and the matrix
exponential and logarithm all come from [](gr_mat). The `ca` context
overrides three methods of the generic matrix interface:
multiplication (`ca_mat_mul`), determinants (`ca_mat_det`) and the pivot
search (`ca_mat_find_pivot`), so that the generic algorithms profit from
the specialisations below. A failed zero test makes a generic function
return `GR_UNABLE`, which the `ca_mat` wrappers translate into a return
value 0 or `T_UNKNOWN`.

## Choosing pivots {#pivots}

`ca_mat_find_pivot` searches column $j$ for a nonzero entry in two
passes. The first pass ignores entries that are rational zeros and takes
the remaining entry of least complexity of representation (according to
`ca_cmp_repr`): a simple entry, such as a small integer or an element of
a small number field, is cheap to test and likely to be nonzero, so a
single zero test usually suffices. If that entry is not proved nonzero,
the second pass tests every entry and takes the simplest one that is
proved nonzero. Entries that are proved zero along the way are replaced
by an exact zero, which simplifies later steps. If no entry is proved
nonzero and not all are proved zero, the search returns `T_UNKNOWN` and
the elimination fails: the rank is then unknown.

## Fraction-free elimination {#fflu}

Gaussian elimination with divisions produces fractions whose numerators
and denominators grow. Bareiss's algorithm [@Bar1968] avoids this: with
$a^{(0)} = A$ and $a_{-1,-1}^{(-1)} = 1$ it updates, for $i, j > k$,

$$
a_{ij}^{(k+1)} = \frac{a_{kk}^{(k)} a_{ij}^{(k)} - a_{ik}^{(k)} a_{kj}^{(k)}}{a_{k-1,k-1}^{(k-1)}},
$$

where by Sylvester's identity every $a_{ij}^{(k)}$ is a $(k+1) \times
(k+1)$ minor of the input, so the division is exact (see
[](algorithms/linear-algebra/fraction-free) for the integer case). `ca_mat_fflu` implements
this with `ca_mat_find_pivot` for pivoting; the last pivot is the
determinant (up to sign) and is returned as the denominator. Over a
field such as $\mathbb{Q}(\pi)$ the entries are rational functions in
the generators, the exact divisions keep them polynomial in the input
entries, and the representations stay smaller than with ordinary
elimination. The division is still a division in [](ca), so fractions
can appear when the formal representation does not see that it is
exact.

## Algorithms without zero tests {#division-free}

Some quantities can be computed without any zero test, which makes them
always succeed:

- The Berkowitz algorithm [@Ber1984] computes the characteristic
  polynomial with $O(n^4)$ ring operations and no divisions, from
  products of Toeplitz matrices built from the leading principal
  submatrices; `ca_mat_charpoly` uses it by default and
  `ca_mat_det_berkowitz` takes the constant coefficient.
- Cofactor expansion (`ca_mat_det_cofactor`, `ca_mat_adjugate_cofactor`)
  is used for matrices of size up to 4 or 5.
- The adjugate is obtained from the characteristic polynomial $p(x) =
  \det(xI - A) = x q(x) + p(0)$ as $\operatorname{adj}(A) = (-1)^{n+1}
  q(A)$, by the Cayley–Hamilton theorem (`ca_mat_adjugate_charpoly`).
  Then $A^{-1} = \operatorname{adj}(A)/\det A$ and $A^{-1} B$ need a
  single zero test, of the determinant, which is how `ca_mat_inv` and
  `ca_mat_nonsingular_solve_adjugate` proceed for small matrices or
  matrices with entries in several fields.

The Danilevsky algorithm [@Dan1937], which reduces $A$ to companion form
by similarity transformations with $O(n^3)$ operations, is faster but
needs pivots and therefore zero tests (`ca_mat_charpoly_danilevsky`).
See [](algorithms/linear-algebra/charpoly) for these algorithms over
general rings.

## Rational and number field matrices {#specialised}

Matrices whose entries are all represented as rational numbers are
recognised with `_ca_mat_same_field` and the computation is delegated
to [](fmpz_mat) or [](fmpq_mat): `ca_mat_mul` (for inner dimension at
least 3) and `ca_mat_det` (for size at least 3), with their multimodular
algorithms (see [](algorithms/linear-algebra/multimodular)). Matrices
whose entries lie in one number field $\mathbb{Q}(a)$ are multiplied by
`ca_mat_mul_same_nf` (when the dimensions are at least 3 and the inner
dimension at least 4), which writes the entries as
polynomials in $a$ over a common denominator, multiplies the resulting
integer polynomial matrices with `fmpz_poly_mat_mul` and reduces modulo
the minimal polynomial; for such matrices, determinants and inverses use
LU decomposition, since zero tests in number fields never fail.

## Eigenvalues and the Jordan form {#jordan}

`ca_mat_eigenvalues` computes the characteristic polynomial and its
roots with multiplicities (`ca_poly_roots`, see [](ca_poly/roots)); it
fails if a factor of degree 4 or more with non-rational coefficients
occurs. For an integer or rational matrix it always succeeds, with
algebraic eigenvalues.

Given the eigenvalues $\lambda$ with algebraic multiplicities $m$,
`ca_mat_jordan_blocks` (through `gr_mat_jordan_blocks`) computes the
ranks $r_k = \operatorname{rank}(A - \lambda I)^k$ for $k = 0, 1, \ldots$
until they stabilise at $n - m$. The number of Jordan blocks of size at
least $k$ for $\lambda$ is $r_{k-1} - r_k$, so the sizes of the blocks
form the conjugate of the partition $(r_0 - r_1, r_1 - r_2, \ldots)$
(the transpose of its Ferrers diagram). `ca_mat_jordan_transformation`
then constructs a basis of generalised eigenvectors from nullspaces, and
`ca_mat_diagonalization` checks whether every geometric multiplicity
$\dim \ker(A - \lambda I)$ equals the algebraic multiplicity.

Matrix functions are computed from the Jordan decomposition $A = P J
P^{-1}$: for a Jordan block $J_\lambda$ of size $s$,

$$
f(J_\lambda) = \sum_{k=0}^{s-1} \frac{f^{(k)}(\lambda)}{k!} N^k,
$$

where $N$ is the nilpotent part, so that $f(A) = P f(J) P^{-1}$ only
needs $f(\lambda)$ and its derivatives at the eigenvalues
(`gr_mat_func_jordan`, used by `ca_mat_exp` with $f = \exp$ and by
`ca_mat_log` with the principal logarithm, which fails for singular
matrices). Polynomials in a matrix
are evaluated with the Paterson–Stockmeyer algorithm [@Paterson1973]
(`ca_mat_ca_poly_evaluate`), with about $2\sqrt{\deg f}$ matrix
multiplications.

## Use in FLINT {#use}

- `ca_mat_find_pivot`, `ca_mat_fflu`, `ca_mat_rref_fflu` and
  `ca_mat_nonsingular_solve_fflu` are specific to [](ca_mat); `ca_mat_lu`,
  `ca_mat_rank`, `ca_mat_rref`, `ca_mat_right_kernel` and
  `ca_mat_nonsingular_solve_lu` wrap the generic [](gr_mat) algorithms.
- `ca_mat_charpoly`, `ca_mat_det`, `ca_mat_adjugate` and `ca_mat_inv`
  choose between the methods above; `ca_mat_mul` and
  `ca_mat_mul_same_nf` specialise multiplication.
- `ca_mat_eigenvalues`, `ca_mat_diagonalization`, the Jordan form
  functions, `ca_mat_exp` and `ca_mat_log` build on root finding and
  nullspaces.
- `qqbar_eigenvalues_fmpz_mat` computes eigenvalues of integer matrices
  directly as algebraic numbers.
