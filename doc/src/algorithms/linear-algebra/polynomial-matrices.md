---
title: Matrices of polynomials
summary: How FLINT multiplies matrices over $\Z[x]$ and $(\Z/n\Z)[x]$ — classical products of polynomials, Kronecker substitution into one integer matrix product, evaluation and interpolation, and transforming each entry once with the small-prime FFT — and how it computes determinants, ranks, inverses and solutions by fraction-free elimination or by evaluation and interpolation.
complexity: For $d \times d$ matrices of polynomials of degree less than $\ell$, multiplication costs $O(d^\omega \ell)$ products over the coefficient ring plus $O(d^2 \ell \log \ell)$ for evaluation and interpolation; fraction-free elimination costs $O(d^3)$ operations on polynomials of degree up to $d\ell$.
order: 70
---

A $d \times d$ matrix whose entries are polynomials of degree less than
$\ell$ can be seen as a polynomial of degree less than $\ell$ with matrix
coefficients. The two views suggest the two families of fast algorithms:
treat the matrix as a big integer matrix (Kronecker substitution) or as a
polynomial to be multiplied by evaluation and interpolation, which
reduces everything to scalar matrix products. FLINT implements both for
`fmpz_poly_mat` ($\Z[x]$) and `nmod_poly_mat` ($(\Z/n\Z)[x]$), and uses
fraction-free elimination ([](algorithms/linear-algebra/fraction-free))
for determinants, ranks and solving. Asymptotically faster algorithms for
the latter, based on minimal approximant bases [@BL1994; @GJV2003], are
not implemented. Matrices over `gr_poly` rings in `gr_mat` use the generic
algorithms of [](gr_mat) instead.

## Multiplication {#mul}

### Classical {#classical}

`fmpz_poly_mat_mul_classical` and `nmod_poly_mat_mul_classical` compute
each entry as a sum of $k$ polynomial products, $mkn$ polynomial
multiplications in all, with the fast polynomial multiplication of
[](fmpz_poly) and [](nmod_poly). For $2 \times 2$ squares,
`nmod_poly_mat_sqr_classical` uses the formulas
$(A^2)_{00} = a^2 + bc$, $(A^2)_{01} = b(a + d)$ and so on, with five
products.

### Kronecker substitution {#ks}

Substituting $x = 2^s$ for a large enough $s$ turns each polynomial into
an integer from which its coefficients can be read off [@Kro1882], and
turns the matrix product into one integer matrix product.

::: lemma {#lem-ks}
Let $A$, $B$ have polynomial entries with coefficients bounded by $2^a$
and $2^b$ in absolute value and lengths at most $\ell_A$, $\ell_B$, and let
$k$ be the inner dimension. The coefficients of the entries of $AB$ are
bounded by $k \min(\ell_A, \ell_B) 2^{a+b}$. If $s$ exceeds the bit size of
this bound by one (for signed coefficients), the entries of $AB$ are
recovered from those of $A(2^s) B(2^s)$.
:::

`fmpz_poly_mat_mul_KS` packs each entry with `fmpz_poly_bit_pack` at
$s = a + b + [\text{signed}] + \lceil \log_2 \min(\ell_A, \ell_B) \rceil
+ \lceil \log_2 k \rceil$ bits, multiplies the integer matrices with
`fmpz_mat_mul`, and unpacks. `nmod_poly_mat_mul_KS` does the same with
coefficients in $[0, n)$, so $s = 2 \log_2 n + \log_2 \min(\ell_A, \ell_B)
+ \log_2 k$ (rounded up), and reduces the unpacked coefficients modulo
$n$. All the multimodular and floating-point machinery of `fmpz_mat_mul`
([](algorithms/linear-algebra/integer-matrix-multiplication)) then
applies to the packed matrices, whose entries have about $\ell s$ bits.

### Evaluation and interpolation {#interpolation}

Over a field with at least $\ell_A + \ell_B - 1$ elements, the product is
determined by its values at that many points. `nmod_poly_mat_mul_interpolate`
evaluates every entry of $A$ and $B$ at the points $0, 1, \ldots,
\ell_A + \ell_B - 2$ with a subproduct tree
([](algorithms/polynomials/subproduct-tree)), multiplies the
$\ell_A + \ell_B - 1$ pairs of scalar matrices with `nmod_mat_mul`, and
interpolates each entry of the product. This costs
$O(\ell \, d^\omega)$ operations modulo $n$ for the products plus
$O(d^2 M(\ell) \log \ell)$ for evaluation and interpolation, so the fast
matrix products of `nmod_mat` dominate for large $d$. It requires $n$ to
be a prime at least $\ell_A + \ell_B - 1$; it throws if $n$ is
smaller.

### Transformed entries {#fft-small}

`nmod_poly_mat_mulmid_fft_small` computes a window of coefficients
$[z_l, z_h)$ of each entry of $AB$ by transforming every entry of $A$ and
$B$ once in the transformed polynomial ring of `gr`
(`gr_ctx_init_gr_poly_transformed_repr`), built on the small-prime FFT of
[](fft_small). Each entry of the product is a sum of $k$ pointwise
products in the transformed domain, followed by one inverse transform
and Chinese remaindering restricted to the window. The transforms of the
$mk + kn$ inputs are computed once instead of $mkn$ times, so for long
polynomials the cost is close to $(mk + kn + mn)$ transforms plus $mkn$
pointwise products, the analogue of
[](algorithms/linear-algebra/integer-matrix-multiplication#transformed).
The conversions and the output entries are distributed over threads. If
the storage for the transforms exceeds the budget of [](fft_small), the
product is computed in row and column blocks. The function returns 0 if
`fft_small` is unavailable or no transform plan exists for the
accumulation bound.

### Choice {#choice}

`fmpz_poly_mat_mul` uses the classical algorithm if any dimension is
less than 8 and Kronecker substitution otherwise; `fmpz_poly_mat_sqr`
likewise. `nmod_poly_mat_mul`, with $d$ the smallest dimension and
$\ell$ the smaller of the maximal lengths of the entries of $A$ and $B$:

| Condition | Method |
|---|---|
| $d \ge 2$, $\ell \ge 2d$ and $\ell \ge c$ | `nmod_poly_mat_mulmid_fft_small` (if it succeeds) |
| $d < 10$ | classical |
| $n$ prime of more than $\beta/4$ bits, $n \ge \ell_A + \ell_B - 1$, $d > 60 + \sqrt{\ell}$ | evaluation and interpolation |
| $\max(\ell_A, \ell_B) > 128$ | classical |
| otherwise | Kronecker substitution |

Here $c = 64$, or 32 when $n$ has more than 40 bits, doubled for $d = 2$.
`nmod_poly_mat_sqr` follows the same rules without the transformed
method, with the interpolation threshold $d > 80 + \sqrt{\ell}$.

### Truncated products and powers {#powers}

`fmpz_poly_mat_mullow` and `fmpz_poly_mat_sqrlow` multiply and truncate
entrywise; `fmpz_poly_mat_pow` and `nmod_poly_mat_pow` use binary
powering, `fmpz_poly_mat_pow_trunc` binary powering with truncation after
each step, and `fmpz_poly_mat_prod` a balanced product tree of the
factors.

## Determinants {#det}

By Laplace expansion, the determinant of a $d \times d$ matrix with
entries of length at most $\ell$ has degree at most $d(\ell - 1)$, so it
is determined by its values at $d(\ell - 1) + 1$ points.

- For $d \le 2$ an explicit formula is used, and for $3 \le d \le 14$
  fraction-free elimination (`fmpz_poly_mat_det_fflu`,
  `nmod_poly_mat_det_fflu`), whose last pivot is the determinant up to the
  sign of the row permutation; the pivots are chosen of low degree.
- For $d \ge 15$, `fmpz_poly_mat_det_interpolate` evaluates the matrix at
  the integers $-\lfloor L/2 \rfloor, \ldots$, $L = d(\ell - 1) + 1$,
  computes the $L$ integer determinants with `fmpz_mat_det`, and
  interpolates with `fmpz_poly_interpolate_exact`.
  `nmod_poly_mat_det_interpolate` evaluates at $0, 1, \ldots, L - 1$, uses
  `nmod_mat_det` and `nmod_poly_interpolate_nmod_vec`, and falls back to
  fraction-free elimination when $L$ exceeds the modulus. Interpolation
  needs the differences of the points to be invertible, so the modulus
  should be a prime.

The interpolation method costs $L$ scalar determinants,
$O(d^{\omega + 1} \ell)$ operations modulo a prime, against $O(d^3)$
operations on polynomials of degree up to $d\ell$ for fraction-free
elimination.

## Rank, solving and inverses {#solve}

`fmpz_poly_mat_rank`, `fmpz_poly_mat_solve`, `fmpz_poly_mat_inv`,
`fmpz_poly_mat_rref` and `fmpz_poly_mat_nullspace`, and their
`nmod_poly_mat` analogues, are built on fraction-free elimination over the
polynomial ring, with exact polynomial divisions. Solutions and inverses
are returned with a polynomial denominator: `fmpz_poly_mat_solve(X, den,
A, B)` gives $AX = \mathit{den} \cdot B$ with $\mathit{den} = \pm \det A$.
By [](algorithms/linear-algebra/fraction-free#thm-sylvester) all
intermediate entries are minors of the input, of degree at most
$d(\ell - 1)$, so the cost is $O(d^3)$ multiplications and exact
divisions of polynomials of that degree, $O(d^3 M(d \ell))$ coefficient
operations over $\Z/n\Z$. The rank is the number of pivots.

## Use in FLINT {#use}

| Task | $\Z[x]$ | $(\Z/n\Z)[x]$ |
|---|---|---|
| Multiplication | `fmpz_poly_mat_mul`, `fmpz_poly_mat_mul_classical`, `fmpz_poly_mat_mul_KS` | `nmod_poly_mat_mul`, `nmod_poly_mat_mul_classical`, `nmod_poly_mat_mul_KS`, `nmod_poly_mat_mul_interpolate`, `nmod_poly_mat_mulmid_fft_small` |
| Squaring, powers | `fmpz_poly_mat_sqr`, `fmpz_poly_mat_pow`, `fmpz_poly_mat_pow_trunc`, `fmpz_poly_mat_prod` | `nmod_poly_mat_sqr`, `nmod_poly_mat_pow` |
| Determinant | `fmpz_poly_mat_det`, `fmpz_poly_mat_det_fflu`, `fmpz_poly_mat_det_interpolate` | `nmod_poly_mat_det`, `nmod_poly_mat_det_fflu`, `nmod_poly_mat_det_interpolate` |
| Elimination | `fmpz_poly_mat_fflu`, `fmpz_poly_mat_rref`, `fmpz_poly_mat_rank`, `fmpz_poly_mat_nullspace` | `nmod_poly_mat_fflu`, `nmod_poly_mat_rref`, `nmod_poly_mat_rank`, `nmod_poly_mat_nullspace` |
| Solving | `fmpz_poly_mat_solve`, `fmpz_poly_mat_solve_fflu`, `fmpz_poly_mat_inv` | `nmod_poly_mat_solve`, `nmod_poly_mat_solve_fflu`, `nmod_poly_mat_inv` |
