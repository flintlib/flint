---
title: Dot products modulo a word-size integer
summary: How FLINT computes $\sum a_i b_i \bmod n$ with a single final reduction, accumulating the exact sum in one, two or three words, or in SIMD lanes, depending on $n$, the length and the processor.
complexity: One product and one to three word additions per term, plus one reduction of at most three words per dot product.
order: 10
---

Consider a dot product modulo a word-size integer: $n \ge 2$ is the
modulus, $a_i, b_i \in [0, n)$ are reduced entries, $\ell$ is the length, $\beta = \bits$ and

$$
S = \sum_{i=0}^{\ell - 1} a_i b_i
$$

is the unreduced dot product. FLINT never reduces modulo $n$ inside the
loop. It computes $S$ exactly, or a value congruent to $S$ that fits in
fewer words, and reduces once at the end with `NMOD_RED` (one word,
[Barrett](algorithms/arithmetic/barrett#nmod)), `NMOD2_RED2` or `NMOD_RED3`
(two or three words,
[Möller–Granlund](algorithms/arithmetic/moller-granlund)). This is the technique of
*delayed modular reduction*, systematically used in exact linear algebra
since FFLAS [@DGP2002; @DGP2008]. Those libraries compute $S$ in
floating-point arithmetic, which is exact while $S < 2^{53}$. FLINT
works with integer words, so the question is how many words $S$ needs.

## How many words {#limbs}

::: lemma {#lem-bound}
$0 \le S \le \ell (n - 1)^2$, and the upper bound is attained when all
entries equal $n - 1$.
:::

The proof is immediate.

::: theorem Limb count {#thm-limbs}
Let $k \ge 1$. Every dot product of length $\ell$ of reduced entries fits
in $k$ words, that is $S < 2^{k \beta}$, if and only if

$$
n \le \left\lceil \frac{2^{k\beta/2}}{\sqrt{\ell}} \right\rceil.
$$
:::

::: proof
By [](#lem-bound), the condition is $\ell (n - 1)^2 < 2^{k\beta}$, that
is $n - 1 < x$ with $x = 2^{k\beta/2} / \sqrt{\ell}$. If $x$ is an integer,
this says $n \le x = \lceil x \rceil$. Otherwise it says
$n - 1 \le \lfloor x \rfloor$, that is $n \le \lfloor x \rfloor + 1 = \lceil x \rceil$.
:::

Some consequences, for $\beta = 64$:

- Since $\ell < 2^{\beta}$ and $(n - 1)^2 < 2^{2\beta}$, three words always
  suffice. If $n \le 2^{\beta/2}$, two words always suffice.
- For $\ell = 1$, one word suffices exactly when $n \le 2^{32}$.
- The thresholds for one and two words are, for example:

| $\ell$ | one word: $n \le$ | two words: $n \le$ |
|---:|---:|---:|
| 2 | 3037000500 | 13043817825332782213 |
| 8 | 1518500250 | 6521908912666391107 |
| 10 | 1358187914 | 5833372668713515885 |
| 100 | 429496730 | 1844674407370955162 |
| 1000 | 135818792 | 583337266871351589 |
| $10^6$ | 4294968 | 18446744073709552 |

`_nmod_vec_dot_params` contains these thresholds as constants for
$2 \le \ell \le 11$ (for $\beta = 32$ and $\beta = 64$). They agree with
[](#thm-limbs) in every case. For longer lengths it computes the three-word
product $\ell (n - 1)^2$ directly, as does `_nmod_vec_dot_bound_limbs`.

## Choice of method {#methods}

`_nmod_vec_dot_params` returns one of the values of `dot_method_t`. The
values are ordered so that the number of words of $S$ can be read off the
method, which `_nmod_vec_dot_bound_limbs_from_params` does. They are tried
in this order:

1. `_DOT0` if $\ell = 0$ or $n = 1$. The result is 0.
2. `_DOT_POW2` if $n$ is a power of 2 (see [](#one-limb)).
3. If $\ell \le 11$: `_DOT1`, `_DOT2` or `_DOT3` according to the number of
   words, by [](#thm-limbs). These are handled by unrolled inline code.
4. If $n \le 2^{\beta/2}$: `_DOT1` if one word suffices. Otherwise, on
   64-bit machines, `_DOT2_SPLIT` if $n \le 1515531528$, and
   `_DOT2_HALF` in the remaining cases (always on 32-bit machines). On
   AArch64, `_DOT2` is used instead of `_DOT2_HALF`, and also instead of
   `_DOT2_SPLIT` when $\ell \ge 160$.
5. If two words suffice: `_DOT_U52`, `_DOT_SPLIT_LIMBS` or `_DOT_U64` if
   available and the length is large enough (see [](#simd)), and otherwise
   `_DOT2`.
6. If three words are needed: `_DOT3_SPLIT_LIMBS` or `_DOT3_U64` if
   available and the length is large enough, otherwise `_DOT3_ACC` if
   $n \le 6521908912666391107$ ($n \le 1518500250$ on 32-bit machines),
   otherwise `_DOT3`.

The parameters are valid for every length up to $\ell$, because each
method's conditions only get weaker as the length decreases. When a SIMD
method is used with a length below its threshold, `_nmod_vec_dot` and
`_nmod_vec_dot_rev` fall back to the scalar code of the same band (`_DOT2`,
`_DOT3_ACC` or `_DOT3`). `_nmod_vec_dot_ptr` and `NMOD_VEC_DOT` always use
this scalar code instead of the SIMD methods for $n > 2^{32}$.

## One word {#one-limb}

If $\ell (n - 1)^2 < 2^{\beta}$, the products and their sum fit in one word
and are added with ordinary arithmetic (`_DOT1`). On 64-bit machines with
AVX2, `_nmod_vec_dot1` keeps four sums in the lanes of a vector register
and multiplies with `vpmuludq`, which forms the 64-bit product of the low
32 bits of each lane. This is exact because one word only suffices for
$\ell \ge 2$ when $n < 2^{32}$. Each lane holds a partial sum of $S$, so
the lanes cannot overflow.

If $n$ is a power of two, $n$ divides $2^{\beta}$, so reduction modulo
$2^{\beta}$ commutes with reduction modulo $n$. The sum can then wrap
around freely and one word always suffices (`_DOT_POW2`). For
$n \le 2^{\beta/2}$ the vectorised `_nmod_vec_dot1` is used; for larger
powers of two, `_nmod_vec_dot_pow2` uses plain scalar code.

## Two and three words {#scalar}

For $n > 2^{\beta/2}$ the products need two words. `_DOT2` computes each
product with `umul_ppmm` and adds it to a two-word sum with
`add_ssaaaa`. It uses two independent sums for the even and odd terms, so
that consecutive additions with carry do not wait for each other, and adds
them at the end.

`_DOT2_HALF` is for $n \le 2^{\beta/2}$ when $S$ needs two words. The
products fit in one word and are added to a two-word sum.

`_DOT3` adds each two-word product to a three-word sum. `_DOT3_ACC` saves
work by adding blocks of 8 products in two words and only the block sums in
three words.

::: proposition {#prop-acc}
The two-word block sums of `_DOT3_ACC` cannot overflow if and only if
$8 (n - 1)^2 < 2^{2\beta}$, that is $n \le 6521908912666391107 \approx 2^{62.5}$
for $\beta = 64$ and $n \le 1518500250 \approx 2^{30.5}$ for $\beta = 32$.
:::

::: proof
A block sum of 8 products is at most $8 (n - 1)^2$, with equality for
entries $n - 1$. The numerical bounds are [](#thm-limbs) with $\ell = 8$
and $k = 2$.
:::

The final reduction of a three-word sum $t_2 B^2 + t_1 B + t_0$, where
$B = 2^{\beta}$, first reduces $t_2$ modulo $n$ and then performs two
two-by-one divisions (`NMOD_RED3`).

## Splitting at bit 56 {#split}

For $n$ somewhat below $2^{\beta/2}$, the two-word sum of `_DOT2_HALF` can be
avoided altogether. On 64-bit machines, `_DOT2_SPLIT` keeps $S$ as
$h \, 2^{56} + l$ in two words $h$ and $l$, moving the high bits of $l$
into $h$ after every 8 terms. At the end, $2^{56}$ is replaced by its
residue $c = 2^{56} \bmod n$, computed in `_nmod_vec_dot_params` and stored
in `pow2_precomp`, and the single word $c h + l$ is reduced. This needs
only one-word additions, which vectorise well.

::: algorithm Dot product with splitting {#alg-split}
**Input:** reduced $a_i, b_i$ for $0 \le i < \ell$, and
$c = 2^{56} \bmod n$.

**Output:** $S \bmod n$.

1. $h \gets 0$, $l \gets 0$.
2. For each full block of 8 terms: add the 8 products $a_i b_i$ to $l$;
   then $h \gets h + \lfloor l / 2^{56} \rfloor$ and
   $l \gets l \bmod 2^{56}$.
3. Add the remaining $r < 8$ products to $l$.
4. Return $(c h + l) \bmod n$, where $c h + l$ is computed as a two-word
   integer: if its high word is 0 it is reduced with a one-word reduction,
   otherwise with a two-word reduction.
:::

The AVX2 version runs this in four lanes, so that each lane receives 8
products per block of 32 terms. It folds and sums the lanes before the
last $r < 4$ products are added.

::: theorem {#thm-split}
Let $m = n - 1$ and let $t$ be the largest number of products added in
step 3, so $t = 7$ for the scalar code and $t = 3$ for the AVX2 code.
[](#alg-split) computes $S \bmod n$ without overflow if

$$
m \le \sqrt{2^{61} - 2^{53}}
\quad\text{and}\quad
\ell \, m^2 < 2^{56} \, (2^{64} - 1).
$$

The first condition is equivalent to $n \le 1515531528 \approx 2^{30.5}$.
At $n = 1515531528$, the second holds for $\ell \le 578721382862847481
\approx 2^{59}$, so it puts no practical constraint on the length.
:::

::: proof
Throughout the algorithm $h \, 2^{56} + l$ equals the sum of the products
added so far, since step 2 only moves multiples of $2^{56}$ from $l$ to
$h$. At the start of each block $l < 2^{56}$, so after 8 more products
$l \le 2^{56} - 1 + 8 m^2$. This is below $2^{64}$ if and only if
$m^2 \le 2^{61} - 2^{53}$, which gives the first condition. Since
$\lfloor \sqrt{2^{61} - 2^{53}} \rfloor = 1515531527$, it is equivalent to
$n \le 1515531528$. The same argument covers step 3 and each lane of the
AVX2 code, where at most 8 products are added between two folds.

In step 4, $h \le \ell \, m^2 / 2^{56} < 2^{64}$ by the second condition,
so $h$ does not overflow, and $l < 2^{64}$. Since $c \le m$,

$$
c h + l < m \, 2^{64} + 2^{64} = n \, 2^{64},
$$

so the high word of the two-word integer $c h + l$ is less than $n$, as
required by the two-word reduction (when it is 0, a one-word reduction
gives the same result). Finally
$c h + l \equiv h \, 2^{56} + l = S \pmod n$.
:::

The split at 56 bits leaves room for the 8 products of a block, which
determines the bound on $n$; a smaller split point would allow slightly
larger moduli. The table in `src/nmod_vec/dot.c` gives the moduli bounds
for split points from 40 to 63 bits (its length columns refer to an
earlier version that computed $c h + l$ in a single word).

`_nmod_vec_dot_params` selects `_DOT2_SPLIT` whenever $n \le 1515531528$.

### Splitting each product {#half}

On 64-bit machines with AVX2, `_nmod_vec_dot2_half` (but not its `_rev`
and `_ptr` forms) applies a variant of the split to $n \le 2^{32}$. Here a
product can have up to 64 bits, so each product $p$ is itself split:
$\lfloor p / 2^{56} \rfloor$ goes to a high sum $h$ and $p \bmod 2^{56}$ to
a low sum $l$. The low lanes are folded every 256 terms. At the end the
result is again $(c h + l) \bmod n$ with $c = 2^{56} \bmod n$.

::: theorem {#thm-half}
With $m = n - 1$ and $n \le 2^{32}$, `_nmod_vec_dot2_half` with AVX2 is
correct if $\ell \, m^2 < 2^{56} (2^{64} - 1)$. This holds for every
$n \le 2^{32}$ when $\ell \le 72057594105036800 \approx 2^{56}$.
:::

::: proof
The invariant $h \, 2^{56} + l = $ (sum of the products so far) holds as
before. Each of the four lanes of $l$ is below $2^{56}$ after a fold and
receives 64 values below $2^{56}$ per block of 256 terms, so it stays below
$2^{64}$. After the last fold, at most 63 more values go into each lane,
so the horizontal sum of the lanes is below
$4 \cdot 64 \cdot 2^{56} = 2^{64}$. After a final fold and the low parts
of the last three products, $l < 2^{56} + 3 \cdot 2^{56} = 2^{58}$. The
high sum satisfies $h \le \ell \, m^2 / 2^{56} < 2^{64}$, and as in
[](#thm-split) the two-word integer $c h + l$ is below $n \, 2^{64}$, so it
is reduced correctly with a single two-word reduction.
:::

The value of `pow2_precomp` is computed for `_DOT2_HALF` exactly when
FLINT itself was built with AVX2.

## SIMD kernels for moduli above $2^{32}$ {#simd}

For $n > 2^{32}$ the products have more than 64 bits and the vector
multiply instructions of x86-64 cannot form them directly. FLINT splits
the entries into smaller
pieces whose products fit in vector lanes, accumulates the pieces in
separate lanes, and recombines them into a two- or three-word integer
before the final reduction. All three kernels compute $S$ exactly, so they
serve both the two-word and the three-word band. For a general discussion
of modular arithmetic with SIMD instructions, see [@HLQ2016].

### Split limbs (AVX2, AVX-512) {#split-limbs}

`_nmod_vec_dot_split_limbs` uses `vpmuludq` ($32 \times 32 \to 64$ bits)
in 4 lanes (AVX2) or 8 lanes (AVX-512F). Let $b$ be the bit length of
$n - 1$.

- For $b \le 61$, each entry is written $x = x_1 2^s + x_0$ with
  $s = \lceil b/2 \rceil$, so $x_0 < 2^s$ and $x_1 < 2^{b - s}$, and
  $xy = x_0 y_0 + 2^s (x_0 y_1 + x_1 y_0) + 2^{2s} x_1 y_1$. Three
  accumulators receive the three parts. Each part of a term is at most
  $P = \max(l_0^2, 2 l_0 l_1)$, where $l_0 = 2^s - 1$ and
  $l_1 = 2^{b-s} - 1$. After every $F = \min(\lfloor (2^{64} - 2^{32}) / P \rfloor, 2^{40})$
  terms, the high 32 bits of every lane are moved to a second accumulator.
  A lane is below $2^{32}$ after such a fold, so it stays below $2^{64}$
  for the next $F$ terms. The second accumulators are emptied into a
  three-word total every 8192 folds, which keeps their lanes below $2^{45}$.
  For example $F = 2048$ for $b = 52$ and $F = 4$ for $b = 61$.
- For $b \ge 62$, the entries are split into 32-bit halves and the four
  partial products $p < 2^{64}$ of each term are added as $X \gets X + p$
  modulo $2^{64}$ and $X_h \gets X_h + \lfloor p / 2^{32} \rfloor$. The
  sum of the low halves of the $p$ is then $X - 2^{32} X_h \bmod 2^{64}$,
  which is exact as long as a lane receives fewer than $2^{32}$ products.
  The accumulators are emptied every 8192 iterations.

### IFMA kernels (AVX-512) {#ifma}

The AVX512-IFMA instructions `vpmadd52luq` and `vpmadd52huq` multiply the
low 52 bits of each lane and add the low or the high 52 bits of the
104-bit product to an accumulator.

- `_nmod_vec_dot_u52` handles $n \le 2^{52}$, so that whole entries are
  multiplied, as in `nmod_mat_mul_u52`. The low and high halves are kept in
  four independent pairs of accumulators to hide the latency of the
  instructions. Every 126 iterations of 32 terms the pairs are added, and
  the lanes summed into a two-word total $L + 2^{52} H$. At most 508
  halves below $2^{52}$ meet in a lane, and
  $8 \cdot 508 \cdot (2^{52} - 1) < 2^{64}$, so neither the lanes nor
  their sum overflow. The total is $S$, so this kernel requires
  $S < 2^{2\beta}$.
- `_nmod_vec_dot_u64` handles any $n$. It splits the entries into 32-bit
  halves and accumulates the low and high 52-bit halves of the four
  partial products in eight accumulators, emptied into a three-word total
  every 250 iterations:
  $$
  xy = x_0 y_0 + 2^{32} (x_0 y_1 + x_1 y_0) + 2^{64} x_1 y_1 .
  $$

The tails of fewer than 8 terms use masked loads.

### When the kernels are used {#simd-use}

The kernels are compiled only for $\beta = 64$, `split_limbs` with AVX2
and the IFMA kernels with AVX-512F, AVX-512DQ and AVX512-IFMA, and neither
when `FLINT_MACHINE_VECTORS_FORCE_GENERIC` or
`FLINT_MACHINE_VECTORS_STRICT_C` is defined. In the two-word band,
`_nmod_vec_dot_params` tries `_DOT_U52`, then `_DOT_SPLIT_LIMBS`, then
`_DOT_U64`; in the three-word band `_DOT3_SPLIT_LIMBS`, then `_DOT3_U64`.
Each is subject to these parameters:

| Parameter | Defined in | Meaning |
|---|---|---|
| `FLINT_NMOD_VEC_DOT_U52_MIN_LEN` | `flint-mparam.h` | least length for `_DOT_U52` (32 or 40 in the shipped files) |
| `FLINT_NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN` | `flint-mparam.h` | least length for split limbs; 0 disables them (48 on Zen 4 up to 160 on Meteor Lake, 0 on Broadwell, 96 by default) |
| `FLINT_NMOD_VEC_DOT_SPLIT_LIMBS_MAX_BITS` | `flint-mparam.h` | split limbs only for $n \le 2^{\text{this}}$; 64 or more means no limit (60 to 64 in the shipped files; 60 on IFMA machines, where `_DOT_U64` is faster above) |
| `NMOD_VEC_DOT_U64_MIN_LEN` | `nmod_vec.h` | least length for `_DOT_U64`, and only for $n > 2^{52}$ (80) |
| `NMOD_VEC_DOT3_U64_MIN_LEN` | `nmod_vec.h` | least length for `_DOT3_U64` (48) |
| `NMOD_VEC_DOT2_SPLIT_MAX_LEN` | `nmod_vec.h` | AArch64 only: `_DOT2` instead of `_DOT2_SPLIT` from this length (160) |

The `_ptr` forms of the kernels exist but are not
selected, since loading `vec2[i][offset]` needs gathers, which measured
slower than the scalar code.

## Use in FLINT {#use}

Dot products with these parameters are the inner loop of classical
polynomial multiplication (`_nmod_poly_mul_classical`, `mullow`,
`mulmid`), power series inversion and division, classical matrix
multiplication ([](nmod_mat/mul)), triangular solving, LU decomposition
with delayed reduction (`nmod_mat_lu_classical_delayed`) and the
characteristic polynomial in `nmod_mat`, and of the `gr` contexts for
$\Z/n\Z$, whose dot product methods are called by the generic code
through `_gr_vec_dot` ([](gr_vec/dot)).

The same bound governs algorithms that delay reductions without forming dot
products. `nmod_mat_lu_classical_delayed` reads the number of words from
the method returned by `_nmod_vec_dot_params` and then updates rows with
the `_nmod_vec_nored` functions in one, two or three words
([](algorithms/linear-algebra/gaussian-elimination#delayed)).
`_nmod_poly_divrem_basecase` does the same with a bit count.

For moduli of several limbs, `mpn_mod` accumulates dot products in
$2\ell - 1$ to $2\ell + 1$ limbs in the same way
([delayed reduction for multiword moduli](algorithms/arithmetic/multiword-delayed-reduction#dot)).
