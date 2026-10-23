---
title: Fixed-precision floating point with nfloat
summary: How nfloat gets its speed - a shallow fixed-size representation selected at run time, kernels specialised by the number of words, fixed-point accumulation for dot products, and fixed-point or exact integer matrix multiplication with rigorous scaling bounds.
complexity: $O(n^2)$ word operations per multiplication of $n$-word numbers; a dot product of length $\ell$ costs about $\ell$ high products of $n + 1$ words; matrix multiplication reduces to Strassen–Waksman fixed-point or multimodular integer matrix products.
order: 10
---

The `nfloat` module implements floating-point numbers whose precision is
a whole number $n$ of words, $1 \le n \le 66$ (64 to 4224 bits on
64-bit machines). It is an example of a *parametric* type of the
generic rings interface [@Joh2025gr]: the element size is fixed when a
context is created, so that elements need no heap memory and pack
contiguously in vectors and matrices, while a single method table serves
all precisions and each operation switches internally to a kernel for
its number of words. For approximate polynomial root-finding and linear
algebra, the generic algorithms run about twice as fast over `nfloat` as
over `arf_t` or MPFR [@Joh2025gr], because of lower overhead in the
basic operations (no allocation, no precision or rounding-mode
arguments, no checks for special values unless enabled) and faster
vector and matrix kernels.

## Representation and basic operations {#representation}

A normal element is $(-1)^s \cdot m \cdot 2^{e - nw}$ with $w$ the word
size and the $n$-word mantissa $m$ normalised, $2^{nw - 1} \le m < 2^{nw}$.
The sign occupies a whole word and special values (zero, infinities,
NaN) are encoded in the exponent, so that the common case needs no bit
manipulation and a single comparison detects all special values. The
exponent is restricted to $\pm 2^{w-2}$, which leaves room to add two
exponents and a correction without overflowing a word, so that overflow
and underflow are detected by one comparison after each operation.

The basic operations dispatch on $n$:

- **Addition and subtraction.** With $|x| \ge |y|$ and
  $\delta = e_x - e_y$, the mantissa of $y$ is shifted right by $\delta$
  bits and truncated to $n$ words; if $\delta \ge nw$, $y$ is negligible.
  For $n \le 4$ there are separate kernels (`_nfloat_add_1`, ...,
  `_nfloat_sub_4`) that work in registers. After a subtraction, a
  cancellation is renormalised by counting leading zeros. The error is
  less than 2 ulp: at most 1 ulp from truncating $y$ and at most 1 ulp
  from truncating the result.
- **Multiplication.** For one word, a $1 \times 1$ word product; for two
  words, a full $2 \times 2$ word product; from three words, the high
  product `flint_mpn_mulhigh_n`, which omits the low partial products at
  the cost of an error of a few units in an extra word, and is then
  normalised by at most one bit. Since the inputs are normalised, the
  product is at least $2^{2nw - 2}$, so a one-bit shift always suffices.
- **Division and square root** use $2/1$ and $3/2$ word divisions with
  precomputed inverses [@GraMol2010] for one and two words, and
  approximate schoolbook division or Newton iterations with a guard
  word otherwise; the guard word is used to certify the direction of
  rounding when directed rounding is requested.

Directed rounding is implemented by computing the truncated result
together with an inexactness flag (or an error bound in units of a
guard word), and adding one ulp in the appropriate cases.

## Dot products {#dot}

Sums of products are the inner loop of polynomial multiplication and
of Gaussian elimination. Rounding each product and each partial sum to
$n$ words, as a loop of `nfloat_mul` and `nfloat_add` would, costs
normalisations and branches per term. `_nfloat_vec_dot` instead
accumulates in fixed point, as `arb_dot` does [@Joh2019dot].

::: algorithm Fixed-point dot product {#alg-dot}
**Input:** $s$ and $x_i, y_i$ for $0 \le i < \ell$, with $n$-word mantissas.

**Output:** an approximation of $s + \sum_i x_i y_i$.

1. Find the largest exponent $E_0 = \max_i (e_{x_i} + e_{y_i})$ (and
   $e_s$), skipping zero terms, and set $E = E_0 + p$ with
   $p = \lfloor \log_2(\ell + 1) \rfloor + 2$ padding bits.
2. Accumulate in an $(n + 1)$-word fixed-point number $S$ with unit
   $2^{E - (n+1)w}$: for each term, with $d = E - (e_{x_i} + e_{y_i})$,
   add or subtract the high product of the mantissas shifted right by
   $d$ bits, computing only the part above the unit; terms with
   $d \ge (n+1)w$ are skipped.
3. Normalise $S$ and round it to $n$ words.
:::

::: proposition {#prop-dot}
In [](#alg-dot), $|S - (s + \sum_i x_i y_i)| \le (c\ell + 1)\, 2^{E - (n+1)w}$,
where $c$ bounds the error of one truncated shifted high product in
units of the accumulator. The accumulator does not overflow, and the
error is at most about $c\ell$ units in the guard word of the result,
relative to the largest term.
:::

::: proof
Each term satisfies $|x_i y_i| < 2^{e_{x_i} + e_{y_i}} \le 2^{E - p}$,
so $|s| + \sum_i |x_i y_i| < (\ell + 1) 2^{E - p} \le 2^{E - 1}$ and the
partial sums fit in the accumulator with a sign. A skipped term is less
than one unit; a computed term is truncated at the unit with an error of
at most $c$ units (the error of the high product plus the truncation of
the shift), and $s$ is truncated once.
:::

The implementation uses $c = 4$ for $n \le 2$ (where the products are
exact two- or three-word products held in registers) and $c = n + 2$
otherwise, and with directed rounding adds or subtracts the bound
$c\ell + 1$ before the final rounding, which makes the result a valid
bound. A first pass over the exponents costs about 25% compared with
maintaining the exponent on the fly, but keeps the inner loop simple.
The complex dot product accumulates the real and imaginary parts in two
such accumulators.

## Matrix multiplication {#matmul}

`nfloat_mat_mul` chooses, by dimension and precision, between the
classical algorithm on top of the dot product and two reductions to
exact or fixed-point arithmetic.

**Fixed-point multiplication** (`nfloat_mat_mul_fixed`). If the
exponents of the entries of $A$ lie in a range of width $\Delta_A$, and
those of $B$ in a range of width $\Delta_B$, then converting $A$ and $B$
to fixed point with a common exponent each, at precision
$p + \Delta_A + \Delta_B$ plus a few guard bits, loses no relative
precision in any entry. The product is then a fixed-point matrix product
`_nfixed_mat_mul`, which needs no normalisation or exponent handling at
all in its inner loops. It uses the classical algorithm with unrolled
dot products of 2 to 8 words, Waksman's algorithm [@Waksman1970], which
saves about half of the multiplications in exchange for additions and
pays off from about 9 words, or Strassen's algorithm [@Str1969] above
dimension 26 to 57. Since fixed-point numbers have no overflow handling,
the inputs are scaled by $2^{-k}$ with $k$ increased until
`_nfixed_mat_mul_bound` certifies that every intermediate value of the
chosen algorithm stays below 1; the same function bounds the error in
ulps. The method fails (and the caller falls back) if the exponent
spread is too large (more than 64 extra bits below 768 bits of
precision, a quarter of the precision above) or if the matrices contain
zeros.

**Block multiplication** (`nfloat_mat_mul_block`). For large matrices
or entries with a wide range of exponents, the inner dimension is split
into blocks within which each row of $A$ and each column of $B$ has
entries of similar magnitude. Each block is scaled row by row and
column by column to an integer matrix and multiplied exactly with
`fmpz_mat_mul`, which uses multimodular or FFT methods; the only errors
are from truncating the entries to integers and the final rounding
[@Joh2019dot]. This turns floating-point matrix multiplication into
integer matrix multiplication with asymptotically fast algorithms.

Complex matrices are multiplied with three real products by the
Karatsuba formula $(AC - BD) + ((A + B)(C + D) - AC - BD)i$ in the
fixed-point and block methods, or by separating real and imaginary parts
into real matrices (`nfloat_complex_mat_mul_reorder`), which is used for
matrices with real entries and for parts with very different scales.
Linear algebra (`nfloat_mat_lu`, triangular solving) uses the recursive
generic algorithms above tuned dimensions, so that most of its work runs
on this matrix multiplication.

## Use in FLINT {#use}

The algorithms described above are implemented in `src/nfloat`: the
arithmetic kernels in `nfloat.c` and `div.c`, the dot products in
`dot.c` (`_nfloat_vec_dot`, `_nfloat_complex_vec_dot`), the fixed-point
arithmetic in `nfixed.c` and the matrix multiplication in `mat_mul.c`.
The elementary functions use the fixed-point kernels of the `mp_real`
module. The `nfloat` contexts are the default choice for fast
approximate arithmetic in the generic modules, for example for
polynomial root-finding with `gr_poly` and approximate linear algebra
with `gr_mat`.

::: functions nfloat_add nfloat_mul _nfloat_vec_dot nfloat_mat_mul nfloat_mat_mul_fixed nfloat_mat_mul_block _nfixed_mat_mul :::
