---
title: SIMD kernels of the element-wise operations
---

On x86-64 with AVX-512 (F, DQ and VL) or AVX2, and on ARM64 with NEON,
`_nmod_vec_reduce`, `_nmod_vec_mul`, `_nmod_vec_scalar_mul_nmod` and
`_nmod_vec_scalar_addmul_nmod` run vector kernels chosen from the size of
$n$. The
instruction set is the one FLINT is compiled for (for instance
`-march=native`, which `configure` sets by default): there is no run-time
detection, and `FLINT_MACHINE_VECTORS_FORCE_GENERIC` or
`FLINT_MACHINE_VECTORS_STRICT_C` turn the kernels off. Below 8 entries the
scalar code runs: inside FLINT functions, where the length varies from call
to call, the kernels lost on such short vectors.

With AVX-512, the last entries are handled by one masked vector rather
than a scalar tail. From 128 entries, when the output is not aligned to 64
bytes and every vector has the same offset modulo 64 bytes, one masked
vector first handles the entries up to the next 64-byte boundary, so that
the main loop reads and writes whole cache lines. On ARM64, from 16
entries, the NEON kernels are interleaved with the scalar code, 4 entries
of each per block of 8 (the scalar multiplier runs beside the NEON one);
the NEON kernels alone handle the rest.

Below, $n$ has $\nu$ bits, $2^{\nu - 1} \le n < 2^{\nu}$, and
$m = \lfloor (2^{\nu + 63} - 1)/n \rfloor$ is the constant of `NMOD_RED`
([](algorithms/arithmetic/barrett#nmod)), on 64-bit machines.

## Reduction of words {#reduce}

`_nmod_vec_reduce` reduces any 64-bit words. Its kernels use products of
32-bit halves only, with the constant
$m_1 = \lfloor m/2^{32} \rfloor = \lfloor (2^{\nu + 31} - 1)/n \rfloor < 2^{32}$
([](algorithms/arithmetic/barrett#nmod)).

| $\nu$ | Kernel | Method | Corrections |
|---|---|---|---|
| 64 | fullword | $a < 2n$ | one |
| 33 to 63 | top32 | $q = \lfloor \lfloor a/2^{32} \rfloor m_1 / 2^{\nu - 1} \rfloor$, then $a - qn$ | one (two at 33) |
| 6 to 32 | fold29 | $t = \lfloor a/2^{\nu + 29} \rfloor (2^{\nu + 29} \bmod n) + (a \bmod 2^{\nu + 29})$, then $q = \lfloor \lfloor t/2^{\nu - 2} \rfloor m_1 / 2^{33} \rfloor$ and $t - qn$ | one |
| 2 to 5 | fold32 | $t = \lfloor a/2^{32} \rfloor (2^{32} \bmod n) + (a \bmod 2^{32})$, then $q = \lfloor s \lfloor m/2^{31} \rfloor / 2^{32} \rfloor$ with $s = \lfloor t/2^{\nu} \rfloor$, and $t - qn$ | two ($2n$, then $n$) |

With 32-bit words, $a$ is a two-word integer and, for $\nu \le 32$, $n$ is
a one-word modulus: these kernels are reductions of
[](algorithms/arithmetic/barrett) at that word size. The top32 kernel is
the estimate of `NMOD_RED` from the top halves of $a$ and $m$; the fold29
kernel is the step of `NMOD_RED2_BOUNDED` after a fold; the fold32 kernel
folds as `n_lll_mod_barrett` does, where `NMOD2_RED2` would reduce the
high half and take a Möller–Granlund step. AVX-512 computes $qn$ for
$\nu \ge 33$ with one `vpmullq`, AVX2 and NEON with two products of 32-bit
halves ($q < 2^{32}$). NEON has no kernel for $\nu \le 5$, where the scalar
code was faster.

::: lemma {#lem-reduce-top32}
Let $33 \le \nu \le 63$, $0 \le a < 2^{64}$ and
$q = \lfloor \lfloor a/2^{32} \rfloor m_1 / 2^{\nu - 1} \rfloor$. Then
$0 \le a - qn < 2n$ for $\nu \ge 34$, and $< 3n$ for $\nu = 33$.
:::

::: proof
Write $a = a_1 2^{32} + a_0$ with $a_0, a_1 < 2^{32}$. From
$(2^{\nu + 31} - n)/n \le m_1 \le 2^{\nu + 31}/n$,
$$
\frac{a}{n} - \frac{2^{32}}{n} - \frac{a_1}{2^{\nu - 1}}
< \frac{a_1 m_1}{2^{\nu - 1}} \le \frac{a_1 2^{32}}{n} \le \frac{a}{n},
$$
so $q \le a/n$, and $q > a/n - 1 - 2^{32}/n - a_1/2^{\nu - 1}$. Both
$2^{32}/n$ and $a_1/2^{\nu - 1}$ are below $2^{33 - \nu}$, so
$a/n - q < 1 + 2^{34 - \nu}$: at most 2 for $\nu \ge 34$ and 3 for
$\nu = 33$, and $a - qn$ is below $2n$, respectively $3n$.
:::

::: lemma {#lem-reduce-fold29}
Let $6 \le \nu \le 32$, $t$ and $q$ as in the table. Then
$t \equiv a \pmod n$, $t < 2^{\nu + 30}$ and $0 \le t - qn < 2n$.
:::

::: proof
The high part $\lfloor a/2^{\nu + 29} \rfloor < 2^{35 - \nu}$ times
$2^{\nu + 29} \bmod n < 2^{\nu}$ is below $2^{35} \le 2^{\nu + 29}$, and the
low part is below $2^{\nu + 29}$, so $t < 2^{\nu + 30}$. This is the
hypothesis of [](algorithms/arithmetic/barrett#thm-bounded) with 32-bit
words, whose constant $\lfloor (2^{\nu + 31} - 1)/n \rfloor$ is $m_1$ and
whose estimate is $q$. That theorem assumes $\nu \le 31$, so that
$t - qn$ fits a 32-bit word; its proof of $0 \le t - qn < 2n$ only uses
$m_1 < 2^{32}$, which also holds for $\nu = 32$, and the kernels compute
$t - qn$ in 64 bits.
:::

For $2 \le \nu \le 5$, the fold is that of
[](algorithms/arithmetic/barrett#two-words) with $2^{32}$ in place of the
word base: $t < 2^{32} n$, which exceeds the bound of
[](algorithms/arithmetic/barrett#thm-bounded). So $s = \lfloor t/2^{\nu} \rfloor < 2^{32}$,
and the argument of [](#lem-reduce-top32) with
$\lfloor m/2^{31} \rfloor = \lfloor (2^{\nu + 32} - 1)/n \rfloor$ gives
$0 \le t - qn < 4n$, hence the two corrections. This constant exceeds
$2^{32}$, so the kernels compute $s \lfloor m/2^{31} \rfloor$ as
$s 2^{32} + s (\lfloor m/2^{31} \rfloor - 2^{32})$. Every product in these
kernels has both factors below $2^{32}$, except $qn$ for $\nu \ge 33$,
which AVX2 and NEON split into two such products.

The full Barrett step of `NMOD_RED`, with the high half of a $64 \times 64$
product emulated from 32-bit products, was 1.3 to 1.6 times slower than
these kernels.

## Products of two vectors {#mul}

`_nmod_vec_mul` multiplies reduced entries. Its kernels:

| $\nu$ | AVX-512 | AVX2 | NEON |
|---|---|---|---|
| 2 to 32 | halflimb | halflimb | halflimb, hybrid from 16 entries |
| 33 to 50 | ifma50 (with AVX512-IFMA), else fp50 | fp50 | fp50 |
| 51 to 62 | red2bounded | scalar code | scalar code |
| 63, 64 | scalar code | scalar code | scalar code |

The halflimb kernels, for $\nu \le 32$, take a 32-bit Barrett step on the
product $p = ab < 2^{2\nu}$, which fits a word:
$q = \lfloor \lfloor p/2^s \rfloor m_1 / 2^{\nu + 31 - s} \rfloor$ with
$s = \nu - 2$, or $s = 30$ for $\nu = 31$ and $s = 32$ for $\nu = 32$, so
that $\lfloor p/2^s \rfloor < 2^{32}$. For $\nu \le 30$, $p < 2^{\nu + 30}$
and [](algorithms/arithmetic/barrett#thm-bounded) with 32-bit words, as
in [](#lem-reduce-fold29), gives $p - qn < 2n$. For $\nu = 31, 32$ the same
argument gives $p/n - q < 1 + 2^s/n + n^2/2^{\nu + 31}$, where the last two
terms, convex in $n$, stay below 3 on $[2^{\nu - 1}, 2^{\nu})$: so
$p - qn < 4n$, and the kernels correct by $2n$, then by $n$.

The fp50 kernels compute $ab - qn$ in $(-n, n)$ with $q$ the rounded
quotient, as `vec4d_mulmod`
([](algorithms/simd/fp-modular-arithmetic#multiplication)), and add $n$ to
a negative result. The ifma50 kernel splits $p = p_1 2^{52} + p_0$ with two
52-bit products, takes $s = \lfloor p/2^{\nu - 2} \rfloor < 2^{\nu + 2}$
and $q = \lfloor s \lfloor m/2^{13} \rfloor / 2^{52} \rfloor$ (for $\nu = 50$,
$\lfloor s \lfloor m/2^{12} \rfloor / 2^{53} \rfloor$, so that the constant
fits 52 bits), and computes $p_0 - qn$ modulo $2^{52}$ as
$p_0 + q(2^{52} - n)$, which is exact since $p - qn < 2n < 2^{52}$: for
$\nu = 50$ this is [](algorithms/arithmetic/barrett#thm-bounded) with
52-bit words, and for $\nu \le 49$ its proof still applies, since
$p < 2^{2\nu} \le 2^{\nu + 49}$, half the bound of that theorem, leaves
room for the halved constant.

The red2bounded kernel, for $51 \le \nu \le 62$, forms the full product
and reduces it as `NMOD_RED2_BOUNDED`
([](algorithms/arithmetic/barrett#bounded)), whose estimate
$\tilde q = \lfloor x m / 2^{65} \rfloor$ leaves a remainder below $2n$;
the high word of $xm$ is computed from three 32-bit products and may be up
to 2 too small, so the quotient may be one less and the remainder is below
$3n$: two corrections. Without AVX-512, this kernel, with its 64-bit
products emulated, was slower than the scalar code.

## Products by a scalar {#scalar}

`_nmod_vec_scalar_mul_nmod` and `_nmod_vec_scalar_addmul_nmod`, for
$c \notin \{0, 1, n - 1\}$, use [Shoup's method](algorithms/arithmetic/shoup)
in several radixes, with $\hat c = \lfloor c 2^{64}/n \rfloor$ computed
once; `_nmod_vec_scalar_addmul_nmod` adds the result and corrects once
more.

| $\nu$ | AVX-512 | AVX2 | NEON |
|---|---|---|---|
| 2 to 32 | shoup32 | shoup32 | shoup32 |
| 33 to 50 | shoup52 (with AVX512-IFMA), else fp50 | fp50 | fp50 (`scalar_mul` only) |
| 51 | shoup52 (with AVX512-IFMA), else shoup64lazy | scalar code | scalar code |
| 52 to 62 | shoup64lazy | scalar code | scalar code |
| 63 | shoup64 | scalar code | scalar code |
| 64 | scalar code | scalar code | scalar code |

In radix $R = 2^{32}$ or $2^{52}$, with $\lfloor \hat c / 2^{64 - k} \rfloor
= \lfloor cR/n \rfloor$ for $R = 2^k$, [](algorithms/arithmetic/shoup#lem-quot)
holds for entries $x < R$: the quotient $q = \lfloor x \lfloor cR/n \rfloor / R \rfloor$
is $\lfloor xc/n \rfloor$ or one less, and $xc - qn < 2n$. With 32-bit
products (shoup32) the difference is computed exactly in 64 bits; with
IFMA (shoup52) it is computed modulo $2^{52}$ as $xc + q(2^{52} - n)$,
which requires $2n < 2^{52}$, hence $\nu \le 51$. In shoup64lazy the high
word of $x \hat c$ is computed from three 32-bit products and may be up to
2 too small, so $xc - qn < 4n$, below $2^{64}$ for $\nu \le 62$, and the
kernel corrects by $2n$, then by $n$; in shoup64, for $\nu = 63$, it is
exact.

The fp50 kernels are the doubles of [](#mul), with $c$ in place of the second
entry; on AVX-512 without IFMA, they were faster than shoup64lazy at 33 to 50
bits. On AVX2, Shoup's method above 32 bits, with emulated 64-bit products, and
on NEON, Shoup's method above 32 bits and the doubles for `scalar_addmul`, were
slower than the scalar code.
