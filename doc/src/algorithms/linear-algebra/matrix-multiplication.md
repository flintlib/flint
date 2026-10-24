---
title: Multiplying matrices over Z/nZ
kind: topic
summary: The methods FLINT uses to multiply matrices modulo a word-size integer, from classical multiplication with delayed reduction to Strassen–Winograd, SIMD kernels, floating-point gemm with Chinese remaindering and bit-packed kernels for tiny moduli, and how `nmod_mat_mul` chooses among them.
order: 5
---

Let $A$ be an $m \times k$ and $B$ a $k \times n$ matrix over $\Z/N\Z$
with $1 \le N < 2^{\bits}$, stored with entries in $[0, N)$. (Throughout
this section $N$ is the modulus and $n$ a dimension.) The product $C = AB$ needs
$mkn$ multiplications and additions of residues with the classical
algorithm, and every fast algorithm in `nmod_mat`, from LU decomposition
to the characteristic polynomial, spends most of its time in such
products. The cost of a modular product is dominated not by the
multiplication but by the reduction modulo $N$, so all methods below
reduce as rarely as possible: they compute exact integer sums of many
products, in machine words or in floating point, and reduce once per
entry or once per block. This *delayed reduction* is the principle of
the FFLAS-FFPACK library [@DGP2002; @DGP2008].

## Generic version {#generic}

`gr_mat_mul` multiplies matrices over any ring through the ring method
`GR_METHOD_MAT_MUL`. The `gr` context for $\Z/n\Z$ installs
`nmod_mat_mul`, so generic matrix code over this ring, such as the
recursive LU decomposition of `gr_mat` or modular composition in
`gr_poly`, uses the methods described here; the contexts
`nmod8` and `nmod32`, whose entries are stored in bytes and 32-bit
words, call the byte and 32-bit kernels of `nmod_mat` directly. Other rings install their own products
([](gr_mat/mul)), and the generic default is the classical algorithm
with one dot product (`_gr_vec_dot`) per entry, the analogue of the
classical method below with the ring's dot product in place of delayed
modular reduction. One-level Strassen-Winograd, Waksman and Rosowski
multiplications (`gr_mat_mul_strassen`, `gr_mat_mul_waksman`,
`gr_mat_mul_rosowski`) are available for rings to install above a
cutoff. The role of multiplication in the other generic algorithms is
described in [](algorithms/linear-algebra/overview).

## The methods {#methods}

| Method | Moduli | Arithmetic | Reductions | Entry point |
|---|---|---|---|---|
| [Classical, delayed](algorithms/linear-algebra/dot-products) | any | integer, 1–3 words per sum | one per entry | `nmod_mat_mul_classical` |
| [Classical, packed](#classical) | $k (N - 1)^2 < 2^{\bits/2}$ | several dot products per word | one per entry | `nmod_mat_mul_classical` |
| [Strassen–Winograd](algorithms/linear-algebra/strassen-winograd) | any | 7 half-size products, 15 additions | in the products and additions | `nmod_mat_mul_strassen` |
| [Integer SIMD](algorithms/linear-algebra/modular-matrix-kernels#u32) | $N < 2^{32}$ | $32 \times 32 \to 64$ bit lanes, symmetric residues | folds when a lane could overflow | `nmod_mat_mul_u32` |
| [IFMA](algorithms/linear-algebra/modular-matrix-kernels#u52) | $N \le 2^{52}$ | 52-bit products, low and high halves | one per entry per block of 256 | `nmod_mat_mul_u52` |
| [Two-limb Karatsuba](algorithms/linear-algebra/modular-matrix-kernels#k52) | $N \le 2^{52}$ | three $32 \times 32 \to 64$ products per product | one per entry per block of 256 | `nmod_mat_mul_k52` |
| [Floating-point mulmod](algorithms/linear-algebra/modular-matrix-kernels#fp50) | $N < 2^{50}$ | doubles, each product reduced | every product, plus folds | `nmod_mat_mul_fp50` |
| [Floating-point gemm](algorithms/linear-algebra/modular-matrix-kernels#blas) | any (64-bit) | `sgemm` or `dgemm`, Chinese remaindering if needed | one per entry | `nmod_mat_mul_blas` |
| [Byte kernels](algorithms/linear-algebra/small-modulus-matrix-multiplication) | $N \le 255$ | table lookups in bytes, `sgemm` above 16 | every few steps | `nmod_mat_mul_u8` |
| [Bit-packed](algorithms/linear-algebra/small-modulus-matrix-multiplication#packed) | $N = 2, 3$ | bitwise operations on 64 entries at once | none | `nmod_mat_mul_u8` |

All of them return the exact product modulo $N$; they differ in speed
only. Every method except the classical and Strassen–Winograd algorithms
and the byte kernels requires 64-bit words.

## Bounds that make delayed reduction exact {#bounds}

All the methods rest on bounds for the size of a sum of products of
residues. With residues in $[0, N)$, a dot product of length $k$ is at
most $k (N - 1)^2$
([](algorithms/linear-algebra/dot-products#limbs)). With residues in the
symmetric range $|x| \le \lceil N/2 \rceil$, which costs one comparison
per entry to set up, the bound drops to about $k N^2 / 4$, four times
less. The integer SIMD kernels and the floating-point methods use the
symmetric representation and the classical algorithm the nonnegative
one.

::: proposition {#prop-exact}
Let $|a_i|, |b_i| \le h$ be integers and let $M$ be a power of two with
$k h^2 < M$. If the partial sums $\sum_{i < j} a_i b_i$, $j \le k$, are
computed in signed integers modulo $2M$, or in a binary floating-point
format with a significand of at least $\log_2 M$ bits, then every partial
sum, and hence the result, is exact.
:::

::: proof
Every partial sum has absolute value at most $j h^2 \le k h^2 < M$, so it
is representable as a signed integer in $(-M, M)$ and as a floating-point
number with $\log_2 M$ significand bits, since all integers of absolute
value at most $M$ are. Induction on $j$ shows that no operation rounds or
wraps around.
:::

For doubles ($M = 2^{53}$) and symmetric residues this allows
$k \lfloor N/2 \rfloor^2 < 2^{53}$, that is $N$ up to about
$2^{27.5}/\sqrt{k}$, or 22 bits for $k = 1000$. For signed 64-bit lanes
($M = 2^{63}$) it allows $N$ up to about $2^{32.5}/\sqrt{k}$ without any
intermediate reduction. Larger moduli need folds (`nmod_mat_mul_u32`),
wider products (`nmod_mat_mul_u52`, `nmod_mat_mul_k52`), a reduction of
every product (`nmod_mat_mul_fp50`) or several products modulo smaller
numbers (`nmod_mat_mul_blas`).

## Classical multiplication {#classical}

`nmod_mat_mul_classical` computes each entry as a dot product with the
parameters of `_nmod_vec_dot_params` for length $k$, which accumulate
the exact sum in one, two or three words and reduce once
([](algorithms/linear-algebra/dot-products)). For dimensions of at least
20 it first transposes $B$, so that both operands of each dot product are
contiguous.

When a dot product fits in one word and all dimensions exceed 10,
several dot products are computed with one word operation. Let $s$ be
the bit length of $k (N - 1)^2$ and $q = \lfloor \bits/s \rfloor$. The
entries $B_{l, j}, \ldots, B_{l, j + q - 1}$ are packed into the word
$w_l = \sum_{t < q} B_{l, j + t} 2^{ts}$, and then

$$
\sum_{l < k} A_{il} w_l = \sum_{t < q} C_{i, j + t} \, 2^{ts}
$$

computed in one-word arithmetic has the unreduced dot products
$C_{i, j+t}$ as its $s$-bit digits: each is at most $k (N - 1)^2 < 2^s$,
so no digit carries into the next and the sum is below
$2^{qs} \le 2^{\bits}$. This pays off when $q \ge 2$, that is
$k (N - 1)^2 < 2^{\bits/2}$, for example for $N < 2^{11}$ and $k < 1000$
on a 64-bit machine.

## Fast multiplication {#fast}

Strassen's algorithm [@Str1969] in Winograd's form [@Win1971] multiplies
$2 \times 2$ block matrices with 7 block products and 15 block additions
instead of 8 products and 4 additions. Applied recursively, it needs
$O(d^{\log_2 7}) = O(d^{2.81})$ operations for dimension $d$. In FLINT,
`nmod_mat_mul_strassen` performs one level and computes the 7 products
with `nmod_mat_mul`, which recurses while the dimensions stay above the
cutoff; see [](algorithms/linear-algebra/strassen-winograd). With the
fast kernels the crossover is high, 256 to 800 depending on the
processor, because the kernels run close to the peak arithmetic rate
while the 15 additions are limited by memory bandwidth.

## Floating point {#floating-point}

A floating-point matrix product (`gemm`) from a tuned BLAS library, or
FLINT's own kernels ([](algorithms/linear-algebra/gemm)), runs at the
peak arithmetic rate of the processor, and by [](#prop-exact) it
computes products of residues exactly as long as
$k \lfloor N/2 \rfloor^2 < 2^{53}$. For larger moduli
`nmod_mat_mul_blas` multiplies modulo up to 12 smaller numbers and
reconstructs the result by Chinese remaindering, which for $k = 1000$
takes 4 `dgemm` calls for 32-bit moduli, 6 for 52-bit moduli and 7 for
64-bit moduli; see
[](algorithms/linear-algebra/modular-matrix-kernels#blas). FLINT's
integer and floating-point SIMD kernels were written to avoid this
multimodular step for $N$ up to $2^{52}$.

## Choosing a method {#choice}

`nmod_mat_mul` makes the choice; the order of the tests and the
parameters are listed in [](nmod_mat/mul#dispatch). In summary, on a
64-bit machine with AVX2, AVX-512 or NEON:

- $N \le 255$: the byte and bit kernels of `nmod_mat_mul_u8` from
  dimension 8.
- $N < 2^{32}$: `nmod_mat_mul_u32`, possibly replaced by a single `dgemm`
  for large matrices when it is exact.
- $2^{32} \le N \le 2^{52}$: `nmod_mat_mul_fp50` (x86-64, $N < 2^{50}$)
  or `nmod_mat_mul_k52` (AArch64, and for $2^{50} \le N \le 2^{52}$).
  On processors with AVX512-IFMA, `nmod_mat_mul_u52` replaces them, and
  also `nmod_mat_mul_u32` for $N \ge 2^{30}$.
- $N > 2^{52}$: classical multiplication with three-word dot products,
  Strassen–Winograd from dimension 200, and `nmod_mat_mul_blas` with
  several `dgemm` calls from dimension 450; on x86-64 with AVX2, when $A$
  has at most 4 rows (8 with IFMA), $k \ge 2$ and $B$ at least 16
  columns, the few-rows code of the next item.
- Matrices with very few rows or columns, for which the kernels would pad
  most of a register tile with zeros, go to code that streams the other
  operand once ([](algorithms/linear-algebra/modular-matrix-kernels#thin)).
- With one thread and dimensions above
  `FLINT_NMOD_MAT_MUL_SIMD_STRASSEN_CUTOFF`, Strassen–Winograd is applied
  on top of the kernels.

With several threads, the kernels split $C$ into blocks of rows or
columns, `nmod_mat_mul_blas` relies on a threaded `gemm`, and the
classical algorithm distributes rows; Strassen–Winograd is not used
above the kernels.

For moduli of several limbs, the generic ring `mpn_mod` chooses between
classical multiplication with delayed reduction, Waksman's algorithm and
a multimodular method built on `nmod_mat_mul`; see [](mpn_mod/mat-mul)
and [delayed reduction for multiword moduli](algorithms/arithmetic/multiword-delayed-reduction#mat-mul).

## Matrix-vector products {#matvec}

A product $Ab$ with a vector costs $mk$ operations and reads every entry
of $A$ once, so it is limited by memory bandwidth.
`nmod_mat_mul_nmod_vec` computes it as $m$ dot products.
`nmod_mat_nmod_vec_mul` computes $aB$ by adding multiples of the rows of
$B$ to accumulators held in vector registers, with the same delayed
reductions as the kernels. When $A$ has a few rows, `nmod_mat_mul` uses
the same code for all of them at once, reading $B$ only once.
