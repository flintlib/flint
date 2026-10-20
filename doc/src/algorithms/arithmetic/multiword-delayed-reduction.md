---
title: Delayed reduction for multiword moduli
summary: How `mpn_mod` evaluates sums of products modulo an $\ell$-limb $n$ exactly in $2\ell - 1$ to $2\ell + 1$ limbs and reduces once, in dot products, polynomial multiplication and division, Euclid's algorithm, Gaussian elimination and matrix multiplication, and how this plugs into the generic algorithms.
complexity: A sum of $k$ products costs $k$ products of $\ell$ limbs, $k$ additions of about $2\ell$ limbs and one reduction, instead of $k$ modular multiplications (each about three products of $\ell$ limbs) and $k$ modular additions.
order: 52
---

A modular multiplication with an $\ell$-limb modulus costs about three
$\ell \times \ell$-limb products: one to form the product and two for the
reduction ([](algorithms/arithmetic/multiword-mulmod)). Most of the time
in linear algebra and polynomial arithmetic is spent in sums of products
$\sum_i a_i b_i$. *Delayed reduction* computes such a sum exactly as an
integer, in a few more limbs than a product, and reduces it once. This
replaces $k$ reductions by one and saves about two thirds of the
multiplications, and modular additions become plain multi-limb
additions. It is the multiword counterpart of the delayed reduction for
word-size moduli ([](algorithms/linear-algebra/dot-products)), which
goes back to FFLAS [@DGP2002; @DGP2008].

Throughout, $B = 2^{\bits}$, the modulus $n$ has $\ell$ limbs and $b$
bits, and residues are integers in $[0, n)$.

## Generic algorithms and specialised kernels {#generic}

The generic algorithms of `gr_vec`, `gr_poly` and `gr_mat` (schoolbook
polynomial multiplication, basecase division, Gaussian elimination,
classical matrix multiplication, ...) are written in terms of ring
operations, and reduce after every operation. They become delayed
algorithms over a particular ring when the ring overrides the low-level
methods they call. `mpn_mod` overrides the vector methods (the dot
products `_gr_vec_dot`, `_gr_vec_dot_rev` and `_gr_vec_dot_strided`),
polynomial multiplication and division, LU decomposition, row reduction
and matrix multiplication; the higher-level generic algorithms (Newton
iteration, half-gcd, block recursive LU, triangular solving, nullspaces,
...) are unchanged and inherit the speed of the kernels. This is the
design described in [@Joh2025gr]. The speedups of `mpn_mod` over the
`fmpz`-based `fmpz_mod` ring come from these basecases, from fewer
temporary allocations, and from cheap additions, which make Karatsuba
and Waksman multiplication worthwhile.

The ring `fmpz_mod` uses a simpler form of the same idea: its dot
products (`_fmpz_mod_vec_dot`) and matrix products (`fmpz_mod_mat_mul`)
compute the exact integer result with `fmpz` arithmetic (`_fmpz_vec_dot`,
`fmpz_mat_mul`) and reduce once at the end, but the intermediate values
are heap-allocated multiprecision integers of variable size.

## How many limbs {#limbs}

::: lemma {#lem-bound}
Let $k \ge 1$ and $0 \le a_i, b_i < n$ for $1 \le i \le k$, and let
$0 \le s_0 < n$. Then

$$
s_0 + \sum_{i=1}^{k} a_i b_i < 2^{2b + \lfloor \log_2 k \rfloor + 1}.
$$

So the sum fits in
$\lceil (2b + \lfloor \log_2 k \rfloor + 1)/\bits \rceil$ limbs.
:::

::: proof
The left-hand side is at most $n - 1 + k(n-1)^2 \le (k+1)(n-1)^2$ for
$n \ge 3$ (for $n \le 2$ the claim is obvious). Since
$k + 1 \le 2^{\lfloor \log_2 k \rfloor + 1}$ and $(n-1)^2 < 2^{2b}$, the
claim follows.
:::

Since $b > (\ell - 1)\bits$, at least $2\ell - 1$ limbs are needed;
$2\ell$ limbs suffice when $2b + \lfloor \log_2 k \rfloor + 1 \le 2\ell\bits$,
that is, when the top limb of $n$ leaves about $(\log_2 k)/2$ bits
free, and $2\ell + 1$ limbs always suffice for $k < B$. The code
computes the exact limb count where it matters for speed (dot products,
schoolbook polynomial multiplication, division) and uses $2\ell + 1$
limbs elsewhere. The final reduction of an $s$-limb sum is
`mpn_mod_set_mpn`, which costs about one modular multiplication.

## Dot products {#dot}

`_mpn_mod_vec_dot` and its reversed and strided variants accumulate the
full products in the number of limbs given by [](#lem-bound). For
$\ell \le 4$ they call generated kernels from `mpn_extras`
(`flint_mpn_dot_tab`, `flint_mpn_dot_rev_tab` and
`flint_mpn_dot_strided_tab`, indexed by the operand sizes and the output
size), which keep the accumulator in registers and use the fixed-size
multiplication macros; for $\ell = 2$ there are dedicated forward and
reversed kernels. For larger $\ell$ each product is formed with
`flint_mpn_mul_n` and added to a $(2\ell + 1)$-limb accumulator. The
initial value and the sign are applied after the reduction. Lengths 1
and 2 use `mpn_mod_mul` and `mpn_mod_fmma` instead.

Every generic algorithm that calls `_gr_vec_dot` benefits: classical
matrix multiplication, triangular solving, power series inversion and
division in the basecase, and many more.

## Polynomial multiplication {#poly-mul}

The coefficients of a product of polynomials are dot products, and all
algorithms of `_mpn_mod_poly_mulmid` compute the integer coefficients of
the product of the residue polynomials and reduce each coefficient once:

- the schoolbook algorithm with the limb count of [](#lem-bound) for
  $k = \min(\mathit{len1}, \mathit{len2})$;
- Karatsuba's algorithm on unreduced, signed coefficients, which grow by
  up to two bits per recursion level, with an output size of
  $2b + 2\lfloor \log_2 k \rfloor + 3$ bits (`_flint_mpn_poly_mulmid_slimbs`);
- Kronecker substitution with $2b + \lfloor \log_2 k \rfloor + 1$ bits
  per coefficient;
- the small-prime FFT, which reconstructs the integer coefficients from
  their residues modulo several 50-bit primes.

See [](mpn_mod/poly-mul) for the choice between them.

## Division with remainder {#divrem}

The schoolbook division of $A$ (length $\mathit{lenA}$) by $B$ (length
$\mathit{lenB}$) computes the quotient coefficients from the top. Each
step reduces the current top coefficient $r$ of the remainder, sets
$q_i = r \cdot \mathrm{lc}(B)^{-1}$ and subtracts $q_i x^i B$ from the
remainder. `_mpn_mod_poly_divrem_basecase_preinv1` keeps the remainder
unreduced: it adds the products $c \cdot B_j$ with $c = n - q_i$, the
residue of $-q_i$, so that all intermediate values are nonnegative
integers. Each coefficient receives at most
$Q = \mathit{lenA} - \mathit{lenB} + 1$ products and starts below $n$,
so by [](#lem-bound) with $k = Q$ it fits in
$\lceil (2b + \lfloor \log_2 Q \rfloor + 1)/\bits \rceil$ limbs. A
coefficient is reduced only when it becomes the top coefficient, and the
low coefficients are reduced at the end. Each quotient coefficient costs
$\mathit{lenB} - 1$ plain products and additions and two modular
operations, instead of $\mathit{lenB} - 1$ modular multiplications and
subtractions.

## Euclid's algorithm {#euclid}

In the Euclidean algorithm for polynomial gcds, almost every step divides
a polynomial of length $m + 1$ by one of length $m$, with a two-term
quotient $q_1 x + q_0$. The remainder coefficients are
$a_i - q_1 b_{i-1} - q_0 b_i$, two products per coefficient.
`_mpn_mod_poly_divrem_q1_preinv1` computes them in one of three ways:

- with one fused reduction per coefficient (`flint_mpn_fmmamod_preinvn`),
  that is, two products and one reduction;
- with the matrix method of preconditioned multiplication for $q_0$ and
  $q_1$, which avoids the full reductions
  ([](algorithms/arithmetic/multiword-mulmod#matrix));
- with a Karatsuba-like trick, using precomputed multipliers
  $Q_0 = -q_0$, $Q_1 = -q_1$ and $Q_0 + Q_1$ for the matrix method. For
  odd $i$,

  $$
  Q_1 b_{i-1} + Q_0 b_i = (Q_0 + Q_1)(b_{i-1} + b_i) - Q_0 b_{i-1} - Q_1 b_i,
  $$

  where $Q_0 b_{i-1}$ was computed for the previous coefficient and
  $Q_1 b_i$ is reused for the next one, which needs only $Q_0 b_{i+1}$
  more. Two coefficients thus cost three products instead of four.

The `fmma` form is used for moduli up to 576 bits or short divisors, the
Karatsuba form for long divisors and large moduli.

## Gaussian elimination {#lu}

`mpn_mod_mat_lu_classical_delayed` stores, next to the reduced matrix, an
unreduced copy of every entry in $2\ell + 1$ limbs. For each pivot
$p$ in row $r$ and column $c$ and each row $i$ below, it computes the
multiplier $e = a_{ic} p^{-1}$ (reduced, stored in $L$) and adds
$(n - e) \cdot u_{rj}$ to the unreduced entry $a_{ij}$ for all columns
$j > c$, where $u_{rj}$ is the reduced pivot row. An entry receives at
most $\min(m, n)$ such products before it is reduced, so by
[](#lem-bound) $2\ell + 1$ limbs suffice for any matrix of fewer than
$B$ rows. Entries are reduced lazily: the next pivot column when it is
searched for a pivot, and the pivot row when it is used.

::: proposition Row updates as one long addition {#prop-long-add}
Let $x = \sum_j x_j B^{sj}$ and $y = \sum_j y_j B^{sj}$ with
$0 \le x_j, y_j$ and $x_j + y_j < B^s$ for all $j$. Then
$x + y = \sum_j (x_j + y_j) B^{sj}$, that is, adding the concatenations
adds the slots.
:::

::: proof
No slot sum reaches $B^s$, so there is no carry from one slot to the
next.
:::

Since no unreduced entry ever exceeds $2\ell + 1$ limbs, the code
forms the products $(n - e) u_{rj}$ for the whole row in a temporary row
with the same layout and adds it to the unreduced row with a single
`mpn_add_n` over $(2\ell + 1)(n - c - 1)$ limbs, which is faster than
one call per entry. `mpn_mod_mat_reduce_row` (row reduction for minimal
polynomials) uses the same unreduced layout.

## Matrix multiplication {#mat-mul}

Classical matrix multiplication uses the dot products above. Waksman's
algorithm halves the number of products by forming products of sums
$(a + b)(a' + b')$ of residues, which are $\ell$-limb (or, if the top bit
of $n$ is set, $(\ell + 1)$-limb) integers; the products are accumulated
with signs in $2\ell + 1$ limbs and each entry is reduced once. The
multimodular algorithm reduces the residues modulo word-size primes,
where `nmod_mat_mul` applies word-size delayed reduction, and reduces the
reconstructed integer entries modulo $n$ once. See [](mpn_mod/mat-mul).

## Accumulation in transform space {#transformed}

The transformed polynomial rings over `mpn_mod`
(`gr_ctx_init_gr_poly_transformed_repr`) take the idea one step further:
products of polynomials are accumulated pointwise in the transform
domain modulo the `fft_small` primes, and only the final sum is
transformed back, reconstructed and reduced modulo $n$. The number of
primes is provisioned for sums of up to `terms_bound` products (with a
factor 4 for the sign bias and the reconstruction), the analogue of the
limb count of [](#lem-bound). This is used for the transformed
representation of polynomial moduli in `gr_poly_preinv_set`.

## Use in FLINT {#use}

- `_mpn_mod_vec_dot`, `_mpn_mod_vec_dot_rev`, `_mpn_mod_vec_dot_strided`:
  [](#dot).
- `_mpn_mod_poly_mulmid` and the individual multiplication algorithms:
  [](#poly-mul).
- `_mpn_mod_poly_divrem_basecase_preinv1`: [](#divrem).
- `_mpn_mod_poly_divrem_q1_preinv1`, used by Euclid's algorithm in
  `_mpn_mod_poly_gcd` and `_mpn_mod_poly_xgcd`: [](#euclid).
- `mpn_mod_mat_lu_classical_delayed`, `mpn_mod_mat_reduce_row`: [](#lu).
- `mpn_mod_mat_mul_waksman`, `mpn_mod_mat_mul_multi_mod`: [](#mat-mul).
- `_fmpz_mod_vec_dot`, `fmpz_mod_mat_mul`, `fmpz_mod_mat_mul_fmpz_vec`:
  exact `fmpz` arithmetic and one reduction ([](#generic)).

::: functions _mpn_mod_vec_dot _mpn_mod_poly_divrem_basecase_preinv1 mpn_mod_mat_lu_classical_delayed :::
