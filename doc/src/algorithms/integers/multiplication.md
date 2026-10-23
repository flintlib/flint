---
title: Multiplying integers
kind: topic
summary: How FLINT multiplies integers of each size — its own assembly basecases, Karatsuba, GMP's Toom–Cook algorithms, and FFT multiplication with fft_small or Schönhage–Strassen — with the algorithms, their costs and the dispatch in flint_mpn_mul.
order: 10
---

Multiplication is the operation on which all other integer arithmetic is
built: division, square roots, GCDs, radix conversion and the evaluation
of series all cost a constant or logarithmic factor times a
multiplication. FLINT multiplies $n$-limb integers with a sequence of
algorithms of decreasing asymptotic cost and increasing overhead: the
schoolbook method, Karatsuba's method, the Toom–Cook family, and fast
Fourier transforms.

Throughout, $B = 2^{\bits}$ is the limb base, $M(n)$ the cost of
multiplying two $n$-limb integers, and the operands are nonnegative limb
arrays $x = \sum_{i < m} x_i B^i$ and $y = \sum_{j < n} y_j B^j$ with
$m \ge n$.

## Overview {#overview}

| Algorithm | Cost for $n \times n$ limbs | Where FLINT uses it (x86-64 with ADX, Skylake thresholds) |
|---|---|---|
| Schoolbook | $n^2$ limb products | FLINT's assembly kernels up to 16 limbs; GMP's basecase inside its Toom–Cook code |
| Karatsuba (Toom-2) | $O(n^{\log_2 3}) = O(n^{1.585})$ | `flint_mpn_mul_toom22` up to 230 limbs |
| Toom-3, 4, 6.5, 8.5 | $O(n^{1.465})$, $O(n^{1.404})$, … | GMP's `mpn_mul`, from about 230 to 500 limbs |
| Small-prime FFT | $O(n \log n)$ in practice | `fft_small` from 500 limbs |
| Schönhage–Strassen | $O(n \log n \log \log n)$ | the [`fft`](fft) module in builds without `fft_small` |

Asymptotically faster algorithms exist, culminating in the $O(n \log n)$
bit-complexity algorithm of Harvey and van der Hoeven [@HvdH2021], but
none is competitive at sizes that fit in memory.

## Schoolbook multiplication {#schoolbook}

The schoolbook method computes $xy = \sum_{j} (x y_j) B^j$ as $n$ passes
of a multiply-and-add of $x$ by one limb (`mpn_addmul_1`), $mn$ limb
products in all. Its cost per limb product is small, and for short
operands every other algorithm loses to it. On x86-64 processors with the
ADX and BMI2 extensions, FLINT has a separate assembly function for every
pair of sizes up to $16 \times 16$ limbs (`flint_mpn_mul_func_tab`), and
squaring functions up to 14 limbs. Each is fully unrolled, keeps operands
in registers where possible and uses the two independent carry chains of
the `adcx` and `adox` instructions, and none has the call overhead and
loop control of GMP's general basecase. ARM64 has FLINT assembly kernels
too (up to 15 limbs for the shorter operand). Other targets use portable
C kernels up to 7 limbs and GMP's `mpn_mul_basecase` up to 16 limbs.

::: summary algorithms/integers/basecase-multiplication :::

## Karatsuba multiplication {#karatsuba}

::: theorem Karatsuba's identity [@KO1962] {#thm-karatsuba}
Let $x = x_1 \beta + x_0$ and $y = y_1 \beta + y_0$. Then
$$
xy = x_1 y_1 \beta^2 + \bigl(x_1 y_1 + x_0 y_0 - (x_1 - x_0)(y_1 - y_0)\bigr) \beta + x_0 y_0 .
$$
:::

::: proof
Expanding, $(x_1 - x_0)(y_1 - y_0) = x_1 y_1 - x_1 y_0 - x_0 y_1 + x_0 y_0$,
so the coefficient of $\beta$ equals $x_1 y_0 + x_0 y_1$, which is the
middle coefficient of $(x_1 \beta + x_0)(y_1 \beta + y_0)$.
:::

With $\beta = B^{\lceil n/2 \rceil}$, a product of $n$-limb numbers thus
needs three products of half the size (the difference $(x_1 - x_0)(y_1 -
y_0)$ is formed from absolute values and a sign) and $O(n)$ additions.

::: theorem {#thm-karatsuba-cost}
Applied recursively down to a basecase of fixed size, Karatsuba's method
multiplies $n$-limb integers in $O(n^{\log_2 3})$ limb operations.
:::

::: proof
The cost satisfies $K(n) \le 3 K(\lceil n/2 \rceil) + c n$. For $n = 2^k$
this unrolls to $K(n) \le 3^k K(1) + c n \sum_{i < k} (3/2)^i = O(3^k) =
O(n^{\log_2 3})$, and monotonicity extends the bound to all $n$.
:::

FLINT's `flint_mpn_mul_toom22` is GMP's `mpn_toom22_mul` (which evaluates
at $0$, $-1$ and $\infty$, the form given above) compiled so that its
recursive calls reach FLINT's basecase kernels instead of GMP's. It also
accepts slightly unbalanced operands, $n \le m < 2n$, splitting at
$\lceil m/2 \rceil$. On x86-64 with ADX, `flint_mpn_mul` uses it for
$n \le 230$ limbs when $5n \ge 4m$; on ARM64 up to 483 limbs (a
threshold measured on Apple M1). Without FLINT assembly it is not used,
and GMP's `mpn_mul` takes over after the basecase.

## Toom–Cook multiplication {#toom}

Karatsuba's method is the case $k = 2$ of a general scheme due to Toom
[@Too1963] and Cook [@Coo1966]. Write $x = X(\beta)$ and $y = Y(\beta)$
for polynomials $X$ and $Y$ of degree $k - 1$ whose coefficients are
blocks of about $n/k$ limbs. Then $xy = W(\beta)$ with $W = XY$ of degree
$2k - 2$, and $W$ is determined by its values at $2k - 1$ points.

::: algorithm Toom–Cook $k$-way multiplication {#alg-toom}
**Input:** $x, y$ of $n$ limbs, an integer $k \ge 2$, and $2k - 1$
distinct points $v_1, \ldots, v_{2k-1}$ among small integers, simple
fractions and $\infty$.

**Output:** $xy$.

1. Split $x$ and $y$ into $k$ blocks of $\lceil n/k \rceil$ limbs,
   giving $X$ and $Y$.
2. Evaluate $X(v_i)$ and $Y(v_i)$ for each $i$ (linear combinations of
   the blocks with small coefficients; at $\infty$, the leading block).
3. Compute the $2k - 1$ products $W(v_i) = X(v_i) Y(v_i)$ recursively.
4. Interpolate the coefficients of $W$ from these values by a fixed
   sequence of additions, shifts and exact divisions by small constants.
5. Return $W(\beta)$, a sum of shifted coefficients with carries.
:::

::: theorem {#thm-toom}
[](#alg-toom) is correct, and applied recursively with a fixed $k$ it
costs $O(n^{\log(2k-1)/\log k})$ limb operations.
:::

::: proof
$W = XY$ has degree $2k - 2$, so it is determined by $2k - 1$ values,
and the interpolation in step 4 solves the Vandermonde system at the
chosen points exactly (the points are distinct, and the projective point
$\infty$ gives the leading coefficient). Since $x = X(\beta)$ and $y =
Y(\beta)$, $W(\beta) = xy$. Steps 1, 2, 4 and 5 cost $O(n)$ limb
operations for fixed $k$, because the evaluation and interpolation
matrices have constant entries. The recurrence $T(n) = (2k - 1)
T(\lceil n/k \rceil) + O(n)$ solves to $O(n^{\log_k (2k-1)})$ as in
[](#thm-karatsuba-cost).
:::

Larger $k$ lowers the exponent ($1.465$ for $k = 3$, $1.404$ for $k = 4$)
but increases the linear overhead, so each $k$ is fastest only in a
range of sizes. Bodrato and Zanoni [@BZ2007] found evaluation points and
interpolation sequences of minimal cost, which GMP implements for
balanced splittings (`toom33`, `toom44`, and `toom6h` and `toom8h` with
half-integer $k$) and for unbalanced ones (`toom32`, `toom42`, `toom43`,
`toom53`, `toom54`, `toom63`), in which $x$ and $y$ are split into
different numbers of blocks. FLINT does not reimplement these: between
the end of its Karatsuba range and the FFT threshold it calls GMP's
`mpn_mul`, which chooses among them with its own tuning [@GMP].

## FFT multiplication {#fft}

For large operands the product is computed as a convolution. Splitting
$x$ and $y$ into $N$ coefficients of $b$ bits gives polynomials whose
product has as coefficients the convolution of the two coefficient
sequences, and a fast Fourier transform computes a convolution of length
$N$ in $O(N \log N)$ ring operations. FLINT has two implementations:

- **Small-prime FFT** (`fft_small`). The convolution is computed modulo
  several primes $p < 2^{50}$ with number-theoretic transforms whose
  butterflies run in double-precision floating point with vector FMA
  instructions, and the coefficients are recovered by the Chinese
  remainder theorem. `mpn_mul_default_mpn_ctx` multiplies integers this
  way. A fixed number of primes suffices for all sizes that fit in
  memory, so the cost behaves like $O(n \log n)$.
- **Schönhage–Strassen** (`fft`). The convolution is computed in
  $\Z/(2^K + 1)\Z$, where $2$ is a root of unity, so that the transforms
  need only shifts and additions, and the pointwise products are
  recursive multiplications. This gives $O(n \log n \log \log n)$
  [@SS1971]. GMP's `mpn_mul` also ends in a Schönhage–Strassen FFT
  (`mpn_mul_fft`) for very large operands [@GKZ2007].

[](algorithms/fft/choosing) compares the transforms, and
[](fft_small/usage) lists the thresholds at which each caller switches.

::: summary algorithms/fft/small-prime-fft :::

::: summary algorithms/fft/schonhage-strassen :::

## Squaring {#squaring}

A square $x^2$ needs only about half the limb products of a general
product in the schoolbook method, since $x_i x_j = x_j x_i$, and only one
forward transform instead of two with FFTs. FLINT therefore detects
squaring: `flint_mpn_sqr` has its own assembly basecases (up to 14
limbs on x86-64 with ADX, 9 on ARM64, none elsewhere), then calls GMP's
`mpn_sqr` (with its Toom squaring variants), and
`fft_small` from `FLINT_FFT_SMALL_SQR_THRESHOLD` limbs. `fmpz_mul` calls
it when both arguments are the same object.

## Short and middle products {#short}

Newton iteration for division and square roots needs only parts of
products: the high half (`flint_mpn_mulhigh_n`), the low half
(`flint_mpn_mullow_n`) or a window (`flint_mpn_mulmid`). For small sizes
these use assembly kernels that skip the unneeded limb products; above
those, Mulders' short product [@Mul2000], which combines a full product
with two recursive short products [@HZ2011]; and above
`FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF` limbs (500 on recent x86-64), an
FFT product of only the needed range (`_mpn_ctx_mpn_mul_range`), which
costs about as much as a full product of the output length. A middle
product is the transposed operation of a product, as explained in
[](algorithms/polynomials/middle-product); `_mpn_ctx_mpn_mul_range`
computes a window shorter than the longer operand with transposed
transforms, at the cost of a product of the shorter operand and the
window ([](algorithms/integers/short-products#transposed-fft)).

::: summary algorithms/integers/short-products :::

## Use in FLINT {#use}

`flint_mpn_mul(r, x, m, y, n)` (for $m \ge n \ge 1$, with $r$ not
overlapping the inputs, see [](mpn_extras/mul)) is an inline function. If a basecase kernel
exists for $(m, n)$ it calls it; otherwise `_flint_mpn_mul` chooses, with
the thresholds of `flint-mparam.h` for the target (the values given are
those of x86-64 processors from Skylake on):

1. $n = 1$: `mpn_mul_1`.
2. $n \le 230$ and $5n \ge 4m$: `flint_mpn_mul_toom22`.
3. $n < T$, unless $m + n > 2T$ and $3n > T$, where $T =$
   `FLINT_FFT_SMALL_MUL_THRESHOLD` $= 500$: GMP's `mpn_mul`.
4. Otherwise: `mpn_mul_default_mpn_ctx` of `fft_small`.

The exception in step 3 sends unbalanced products to the FFT earlier,
because GMP's cost grows like $(m/n) M(n)$ while that of the FFT depends
on $m + n$. Without `fft_small`, step 4 is replaced by
`flint_mpn_mul_fft_main` of the `fft` module from
`FLINT_FFT_MUL_THRESHOLD` $= 32000$ limbs. `flint_mpn_mul_n` (balanced
operands) and `flint_mpn_sqr` follow the same pattern, and
`flint_mpn_mul_shared` computes two products $ax$ and $bx$ with a common
factor, sharing the transform of $x$ when both are in the FFT range.

The combined effect, for balanced $n \times n$ products measured on the
machine that built this documentation (an x86-64 server with AVX-512,
shared with other jobs, so the ratios vary by about 10%):

| $n$ (limbs) | 4 | 8–16 | 32 | 100 | 200–400 | 800 | 2000 | 10000–50000 |
|---|---|---|---|---|---|---|---|---|
| time of GMP's `mpn_mul_n` / time of `flint_mpn_mul_n` | 2.2 | 1.2 | 1.5 | 1.2 | 1.0–1.1 | 1.5 | 2.6 | 2.1 |

Every product of large integers in FLINT reaches `flint_mpn_mul`:

- `fmpz_mul`, `fmpz_addmul`, `fmpz_pow_ui` and the other `fmpz`
  functions ([](fmpz/mul)), and through them `fmpq` and the matrix and
  polynomial modules over $\Z$ and $\mathbb{Q}$;
- Kronecker substitution in `fmpz_poly` and `nmod_poly`, which packs
  polynomials into large integers ([](algorithms/polynomials/kronecker));
- FLINT's own division, square root and radix conversion
  ([](algorithms/integers/division),
  [](algorithms/integers/radix-conversion));
- the floating-point and ball arithmetic of `arf`, `arb` and `mp_real`,
  and binary splitting ([](algorithms/integers/binary-splitting)).

::: functions flint_mpn_mul flint_mpn_mul_n flint_mpn_sqr fmpz_mul :::
