---
title: Matrix multiplication kernels modulo n
summary: The blocked SIMD kernels `nmod_mat_mul_u32`, `_u52`, `_k52` and `_fp50`, which multiply matrices modulo moduli up to $2^{52}$ in one pass with delayed reduction, and the floating-point route of `nmod_mat_mul_blas` with Chinese remaindering.
complexity: $mkn$ vector multiply-adds (u32, u52 in single mode), two (u52), three (k52) or about six (fp50) per product, plus $O(mn \lceil k/256 \rceil)$ reductions; `nmod_mat_mul_blas` costs one to twelve floating-point products.
order: 20
---

The blocked kernels described below are those behind `nmod_mat_mul`
for moduli $1 \le N \le 2^{52}$; the floating-point method for larger
moduli is described after them.
Throughout, $A$ is $m \times k$, $B$ is $k \times n$, and the entries are
residues in $[0, N)$. The kernels have the same structure as
floating-point `gemm` kernels ([](algorithms/linear-algebra/gemm)); what
differs is how the products are represented and when they are reduced.
The bounds that make delayed reduction possible are summarised in
[](algorithms/linear-algebra/matrix-multiplication#bounds).

## Common structure {#structure}

The four kernels `nmod_mat_mul_u32`, `nmod_mat_mul_u52`,
`nmod_mat_mul_k52` and `nmod_mat_mul_fp50` are instances of one template,
`src/nmod_mat/mul_blocked_templ.h`, organised as in GotoBLAS and BLIS
[@GvdG2008; @VZvdG2015]:

- $C$ is computed in blocks of $k_c = 256$ steps of the inner dimension,
  $m_c = 96$ rows of $A$ and $n_c = 2048$ columns of $B$. Each block of
  $A$ and $B$ is *packed* into a contiguous buffer, in the order in which
  the micro-kernel reads it, after converting the entries to the kernel's
  representation.
- The micro-kernel keeps an $m_r \times n_r$ tile of $C$ in vector
  registers as $m_r \times (n_r / v)$ accumulators of $v$ lanes. Each of
  the $k_c$ steps broadcasts one entry of $A$ per row and multiplies it
  by $n_r / v$ vectors of $B$.
- At the end of each block of $k_c$ steps, the accumulators are reduced
  to residues in $[0, N)$ and stored in $C$, after adding the previous
  contents of $C$ for all blocks but the first. $C$ itself thus carries
  the partial sums from one block to the next, and no wide temporary
  copy of $C$ is needed.
- With several threads, $C$ is split along its longer side into blocks of
  whole tiles, one per thread, each computed by the serial code; at most
  one thread is used per $5 \cdot 10^5$ products $mkn$. Every thread
  packs the shared operand itself.

The cost is $mkn$ multiply-accumulate steps of the kernel and
$mn \lceil k / 256 \rceil$ reductions. Since a reduction is a dozen
vector operations, it costs a few percent of the total.

## Integer kernel for $N < 2^{32}$ {#u32}

`nmod_mat_mul_u32` lifts each entry $a$ to the signed 32-bit integer
$\hat a \equiv a$ with $|\hat a| \le \lfloor N/2 \rfloor$, and multiplies
with the widening $32 \times 32 \to 64$-bit instructions available on
all SIMD instruction sets (`vpmuldq` on x86-64, `smlal` on AArch64). A
product is then at most $P = \lfloor N/2 \rfloor^2 < 2^{62}$ in absolute
value, and the accumulators are signed 64-bit lanes.

When an accumulator could overflow, it is *folded*: with
$c = 2^{32} \bmod N$, write $x = h 2^{32} + \ell$ with $0 \le \ell < 2^{32}$
(so $h$ is the arithmetic right shift of $x$ by 32) and replace $x$ by
$x' = hc + \ell \equiv x \pmod N$. This costs a shift, a mask, one
widening multiplication and an addition.

::: proposition {#prop-fold}
If $|x| < 2^{63}$, then $|x'| \le R_1 = 2^{31} c + 2^{32}$. A second
fold of $x'$ gives a value of absolute value at most
$R_2 = (\lfloor R_1 / 2^{32} \rfloor + 1) c + 2^{32}$.
:::

::: proof
$|h| \le 2^{31}$ and $0 \le \ell < 2^{32}$, so $|hc + \ell| \le 2^{31} c + 2^{32}$.
For the second fold, $|h'| \le \lfloor R_1 / 2^{32} \rfloor + 1$ by the
same decomposition.
:::

The kernel uses one fold if $R_1 \le 2^{40}$, which happens when $c$ is
small, as for $N = 2^{31} - 1$ or $N$ close to $2^{32}$, and two folds
otherwise; let $R$ be the corresponding bound. An accumulator starts each
run from a fold or from an entry of $C$, so it stays below $2^{63}$ in
absolute value for

$$
\min\left(256, \left\lfloor \frac{2^{63} - 1 - \max(R, N)}{P} \right\rfloor\right)
$$

steps, after which it is folded again. This *cadence* is 256, so that
no fold is needed inside a block, for $N$ up to about $2^{28.5}$; it is
128 at 29 bits, 32 at 30 bits, 8 for $N = 2^{31} - 1$ and 2 for $N$ close
to $2^{32}$. At the end of a block, $x$ is reduced exactly: the quotient
$q = \operatorname{round}(x / N)$ is computed in double precision, and
$x - qN$, which lies in $(-N, 2N)$, is corrected by one addition or
subtraction of $N$. The kernel is correct for every $N < 2^{32}$ with no
case analysis; it slows down only near $2^{32}$, where the folds become
frequent.

The same kernel serves `_nmod_mat_mul_u32` on arrays of 32-bit entries,
such as those of the `nmod32` contexts of `gr`.

## IFMA kernel for $N \le 2^{52}$ {#u52}

The AVX512-IFMA instructions `vpmadd52luq` and `vpmadd52huq` multiply the
low 52 bits of two 64-bit lanes and add the low, respectively high, 52
bits of the 104-bit product to an accumulator. For entries
$a, b < N \le 2^{52}$ this is the exact product $ab = 2^{52} h + \ell$,
with $h, \ell < 2^{52}$. Two accumulators per lane, for the low and high
halves, can absorb $2^{12} - 1 = 4095$ such halves without overflowing 64
bits, so a block of $k_c = 256$ steps needs no reduction at all, and the
entries need no lifting.

For $N \le 2^{26}$ the products are below $2^{52}$, the high halves are
zero, and `nmod_mat_mul_u52` uses a second micro-kernel with only the low
accumulator: one instruction per product and twice as many rows per
tile. In this mode it is faster than `nmod_mat_mul_u32`, which needs a
multiplication and an addition per product; in the two-accumulator mode
it is slower than
`nmod_mat_mul_u32` up to about 30 bits, which is why `nmod_mat_mul` uses
it only from 31 bits on (`FLINT_NMOD_MAT_MUL_U52_MIN_BITS`).

At the end of a block, the two sums $L, H < 2^{61}$ represent
$2^{52} H + L$, which is reduced as follows, with $c = 2^{52} \bmod N$:

1. $h \gets H \bmod N$, by a Barrett step through double precision;
2. $r \gets hc \bmod N$, up to one multiple of $N$, by
   [Shoup multiplication](algorithms/arithmetic/shoup) with the
   precomputed quotient $\lfloor c \, 2^{52} / N \rfloor$, the two
   104-bit products coming from IFMA instructions;
3. $C \gets (r + L) \bmod N$, by a second Barrett step.

The double-precision quotients are off by less than one in each Barrett
step, so one correction on each side suffices.

## Two-limb Karatsuba kernel for $N \le 2^{52}$ {#k52}

Without IFMA, `nmod_mat_mul_k52` reaches 52-bit moduli with 32-bit
multiplications. Each entry is lifted to $|\hat a| \le N/2 < 2^{51}$ and
split as

$$
\hat a = a_1 2^{27} + a_0, \qquad -2^{26} \le a_0 < 2^{26}, \quad |a_1| \le 2^{24},
$$

and a product is
$\hat a \hat b = a_1 b_1 2^{54} + (a_0 b_1 + a_1 b_0) 2^{27} + a_0 b_0$.
The middle coefficient is obtained by Karatsuba's trick, so three
accumulators receive $a_0 b_0$, $(a_0 + a_1)(b_0 + b_1)$ and $a_1 b_1$:
three widening multiply-adds per product instead of four. Each of the
three partial products is at most $(2^{26} + 2^{24})^2 < 2^{52.7}$ in
absolute value, so 256 of them sum to less than $2^{60.7}$ in a signed
64-bit lane: no reduction is needed inside a block.

At the end of a block, each sum $x$ is first brought to $|x| < 0.51 N$:
with $x = 2^{32} h + \ell$, the integers $h$ and $\ell$ are exact as
doubles, and $hc + \ell$ with $c = 2^{32} \bmod N$ is reduced with the
floating-point modular multiplication described next. The three
residues $r_0, r_M, r_2$ are then combined as
$r_0 + (r_M - r_0 - r_2) \cdot (2^{27} \bmod N) + r_2 \cdot (2^{54} \bmod N)$
with two more floating-point modular multiplications, and the result is
reduced to $[0, N)$. The floating-point reduction is what limits the
kernel to $N \le 2^{52}$.

## Floating-point kernel for $N < 2^{50}$ {#fp50}

`nmod_mat_mul_fp50` stores the lifted entries $|\hat a| \le N/2$ as doubles
and reduces every product as it is formed, with the floating-point
modular multiplication of `fft_small` [@HLQ2016]:

$$
p = \hat a \hat b, \qquad q = \operatorname{round}(p \cdot N^{-1}), \qquad
r = (\hat a \hat b - p) - (q N - p),
$$

where $\hat a \hat b - p$ and $qN - p$ are computed exactly with fused
multiply-adds. For $N < 2^{50}$ the result satisfies $r \equiv \hat a \hat b$
and $|r| < \frac98 N$. The values $r$ are summed in a double accumulator,
which is exact while it stays below $2^{53}$. Starting from $|x| \le N$,
this holds for
$\lfloor \frac89 (2^{53}/N - 1) \rfloor$ steps, after which the
accumulator is brought back to $(-0.51 N, 0.51 N)$ by subtracting
$\operatorname{round}(x/N) N$. The cadence is at least 256, so that no
fold is needed inside a block, for $N \le 2^{53}/289 \approx 2^{44.8}$,
and it is 6 for $N$ close to $2^{50}$.

The kernel costs about six floating-point operations per product,
against two integer operations for `nmod_mat_mul_u32` and three for
`nmod_mat_mul_k52`, but it keeps a single accumulator per entry of the
tile where `nmod_mat_mul_k52` needs three. On x86-64 this makes it the
faster of the two up to $2^{50}$; on AArch64, where a widening
multiply-add is a single instruction, `nmod_mat_mul_k52` is faster.

## Floating-point gemm and Chinese remaindering {#blas}

`nmod_mat_mul_blas` uses `flint_sgemm` or `flint_dgemm`, which call an
external BLAS when FLINT is configured with one and FLINT's own kernels
([](algorithms/linear-algebra/gemm)) otherwise. The entries are lifted
to $|\hat a| \le \lfloor N/2 \rfloor$ and converted to floating point,
the product is computed, and each entry of the result is converted back
to an integer and reduced modulo $N$. With $h = k \lfloor N/2 \rfloor^2$:

- if $h < 2^{24}$, one `sgemm` is exact;
- if $h < 2^{53}$, one `dgemm` is exact;
- otherwise the method is multimodular.

In the multimodular case, the entries of $AB$ over the integers, with the
original representatives in $[0, N)$, lie in $[0, k (N-1)^2]$. FLINT
chooses pairwise coprime integers $t_0 > t_1 > \cdots > t_{s-1} \ge 100$,
starting from $t_0 = \lfloor \sqrt{2^{55}/k - 1} \rfloor$ and decreasing,
until $\prod t_i > k (N - 1)^2$. With entries lifted to
$|x| \le \lfloor t_i / 2 \rfloor$, the product modulo $t_i$ is exact in
double precision, since $k \lfloor t_i/2 \rfloor^2 \le k t_i^2 / 4 < 2^{53}$.
At most 12 moduli are allowed; otherwise, for an extremely large $k$, the
function returns 0. Typical counts are 4 moduli for 32-bit $N$, 6 for
52-bit and 7 for 64-bit $N$ at $k = 1000$.

::: algorithm Mixed-radix reconstruction {#alg-crt}
**Input:** residues $u_i = x \bmod t_i$ of an integer
$0 \le x < \prod_i t_i$, for $t_0 < t_1 < \cdots < t_{s-1}$ (sorted
increasingly).

**Output:** $x \bmod N$.

1. Compute the digits $v_0 = u_0$ and, for $i \ge 1$,
   $v_i = (u_i - v_0 - v_1 t_0 - \cdots - v_{i-1} t_0 \cdots t_{i-2})
   (t_0 \cdots t_{i-1})^{-1} \bmod t_i$, using precomputed inverses.
2. Return $\sum_i v_i \, (t_0 \cdots t_{i-1} \bmod N) \bmod N$, accumulated
   in two words and reduced once.
:::

This is Garner's algorithm [@Gar1959], see also [@GG2003, §5.6]. Since
$x = \sum_i v_i t_0 \cdots t_{i-1}$ with $0 \le v_i < t_i$ exactly, the
result is correct. The conversions and the reconstruction are
distributed over the thread pool; before each `gemm`, the worker threads
are returned to the pool so that FLINT's own `gemm` can use them.

`nmod_mat_mul` uses this method for moduli below $2^{32}$ only when one
`dgemm` suffices and the dimension is above
`FLINT_NMOD_MAT_MUL_BLAS_1PASS_CUTOFF`, where a tuned BLAS can beat
`nmod_mat_mul_u32`. For moduli above $2^{52}$, where no SIMD kernel
exists, it uses it from dimension 450 (with one thread), where it beats
classical multiplication combined with Strassen–Winograd.

## Thin matrices {#thin}

The kernels pad the last tile of $C$ to $m_r$ rows (4 to 14) and $n_r$
columns (8 to 24), which wastes most of the work when $A$ has very few
rows or $B$ very few columns. `nmod_mat_mul` therefore uses two other
methods for such shapes when a vector code path exists for $N$.

- *Few rows* (at most 4, or 8 with IFMA): the product is computed as a
  sum of the rows of $B$ multiplied by the entries of $A$, reading each
  row of $B$ once for all rows of $A$. The accumulators stay in vector
  registers (or, for very wide $B$, in the level 1 cache). With IFMA,
  low and high 52-bit halves are accumulated and reduced every 511 rows;
  with AVX2 or AVX-512 and $N \le 2^{32}$, products are split at bit 56
  as in [](algorithms/linear-algebra/dot-products#split) and folded every
  255 rows; otherwise, for $N < 2^{50}$, the floating-point modular
  multiplication is used with a fold every four rows. For the larger
  moduli (with AVX2 or AVX-512; above $2^{52}$ with IFMA) the entries are
  split into two limbs of at most 32 bits whose products are accumulated
  in three to six vector lanes per entry of $C$, as in
  [](algorithms/linear-algebra/dot-products#split-limbs), carried into
  each other every few rows (for $N \le 2^{61}$, every $2^{63-b}$ rows
  or so, $b$ the bit length of $N - 1$) or never (32-bit limbs, larger
  $N$), and each entry is reduced once from up to three words at the end;
  when $B$ takes more than 512 KB (half of a level 2 cache of 1 MB, a
  compile-time parameter), the rows go 16 at a time and the accumulators
  are kept in the level 2 cache in between. `nmod_mat_mul` uses it when $B$
  has at least 6 columns. This is also the algorithm of
  `nmod_mat_nmod_vec_mul` (for $B$ wide enough).
- *Few columns* (2 to 4, or 2 to 8 with IFMA): the columns of $B$ are
  transposed once, and each row of $A$ is multiplied with all of them as
  simultaneous dot products that share the loads of the row, with the
  IFMA dot product kernels of `nmod_vec`
  ([](algorithms/linear-algebra/dot-products#ifma)) or the floating-point
  modular multiplication.
