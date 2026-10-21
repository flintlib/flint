---
title: crt_helpers.h
header: crt_helpers.h
description: Inline fixed-length carry chains and the Chinese remaindering kernels shared by fft_small and the multimodular code.
category: Integers
order: 6
---

This header-only module provides the inline building blocks of Chinese
remaindering for a small, compile-time number of primes: additions and
subtractions of fixed-length limb arrays written as straight-line carry
chains, multiply-accumulate kernels for sums $\sum_j C_j y_j$ of fixed
lengths, the final reduction of such a sum modulo the product of the
primes, and a template that instantiates the complete reconstruction loop
of the [`fft_small`](fft_small) convolution engine. It is included by
`fft_small`, by the multimodular reduction of `mpn_extras`
([](mpn_extras/multi_crt)) for small products of moduli, and by `mpn_mod`,
`mp_real` and `radix` for their `fft_small`-based products.

Only the fixed-length additions are meant for general use. Everything
else is described below as internal: some names (`_mul`, `_madd`) are
short and generic, so the header should only be included in source
files, never in other headers.

## Fixed-length addition and subtraction {#addsub}

::: function
void multi_add_0(ulong FLINT_UNUSED(z[]), const ulong FLINT_UNUSED(a[]))
void multi_add_1(ulong z[], const ulong a[])
void multi_add_2(ulong z[], const ulong a[])
void multi_add_3(ulong z[], const ulong a[])
void multi_add_4(ulong z[], const ulong a[])
void multi_add_5(ulong z[], const ulong a[])
void multi_add_6(ulong z[], const ulong a[])
void multi_add_7(ulong z[], const ulong a[])
void multi_add_8(ulong z[], const ulong a[])

Set $z \gets z + a$ modulo $B^n$ for the $n$-limb arrays $z$ and $a$,
where $n$ is the suffix ($0 \le n \le 8$) and $B = 2^{\bits}$; the carry
out is discarded. `multi_add_0` does nothing.

::: details
On x86-64 with AVX2 and on ARM64 these expand to the multi-limb macros
of [`longlong.h`](longlong) (`add_ssaaaa` and its longer relatives),
which compile to a single chain of add-with-carry instructions;
elsewhere, including 32-bit platforms, the macros are used up to three
limbs and `mpn_add_n` beyond (the macro `CRT_HELPERS_UNROLLED_ADD` tells
which). The arrays must be identical or not overlap.
:::
:::

::: function
void multi_sub_0(ulong FLINT_UNUSED(z[]), const ulong FLINT_UNUSED(a[]))
void multi_sub_1(ulong z[], const ulong a[])
void multi_sub_2(ulong z[], const ulong a[])
void multi_sub_3(ulong z[], const ulong a[])
void multi_sub_4(ulong z[], const ulong a[])
void multi_sub_5(ulong z[], const ulong a[])
void multi_sub_6(ulong z[], const ulong a[])
void multi_sub_7(ulong z[], const ulong a[])
void multi_sub_8(ulong z[], const ulong a[])

Set $z \gets z - a$ modulo $B^n$, the borrow out being discarded.
:::

::: function
void multi_rsub_0(ulong FLINT_UNUSED(z[]), const ulong FLINT_UNUSED(a[]))
void multi_rsub_1(ulong z[], const ulong a[])
void multi_rsub_2(ulong z[], const ulong a[])
void multi_rsub_3(ulong z[], const ulong a[])
void multi_rsub_4(ulong z[], const ulong a[])
void multi_rsub_5(ulong z[], const ulong a[])
void multi_rsub_6(ulong z[], const ulong a[])
void multi_rsub_7(ulong z[], const ulong a[])
void multi_rsub_8(ulong z[], const ulong a[])

Set $z \gets a - z$ modulo $B^n$ (reverse subtraction).
:::

::: macro
CRT_HELPERS_UNROLLED_ADD

1 if the fixed-length additions and subtractions are straight-line
carry chains for all lengths (on x86-64 with AVX2 and on ARM64), and 0
if they call `mpn_add_n` and `mpn_sub_n` beyond three limbs.
:::

## Internal functions {#internal}

::: function
unsigned char _addcarry_ulong(unsigned char cf, ulong x, ulong y, ulong * s)
unsigned char _subborrow_ulong(unsigned char cf, ulong x, ulong y, ulong * s)

Set `*s` to $x + y + c$ (respectively $x - y - c$) modulo $B$ for the
carry (borrow) $c \in \{0, 1\}$ and return the carry (borrow) out. With
AVX2 on 64-bit machines these are the compiler intrinsics
`_addcarry_u64` and `_subborrow_u64`; otherwise they are written with the
two-limb macros of `longlong.h`.
:::

::: function
void _mul(ulong * hi, ulong * lo, ulong y, ulong x)
void _madd(ulong * hi, ulong * lo, ulong y, ulong x)

Set $(\mathit{hi}, \mathit{lo})$ to the two-limb product $xy$, respectively
add $xy$ to the two-limb $(\mathit{hi}, \mathit{lo})$ modulo $B^2$. They
use `unsigned __int128` where available.
:::

### Reconstruction kernels {#kernels}

Through token pasting, the header defines a family of inline functions
for each pair $(n, m)$ with $n \in \{m, m + 1\}$ and $1 \le m \le n \le 8$:

```c
void _big_mul_n_m(ulong r[], ulong t[], ulong C[], ulong y);
void _big_addmul_n_m(ulong r[], ulong t[], ulong C[], ulong y);
void _reduce_big_sum_n(ulong r[], ulong t[], const ulong * limit);
```

(for instance `_big_mul_3_2` and `_reduce_big_sum_3`; the trivial cases
`_big_mul_1_0` and `_big_addmul_1_0` exist as well). Here `C` is an
$m$-limb cofactor and $y$ a limb. `_big_mul_n_m` sets an $n$-limb
accumulator to $Cy$, and `_big_addmul_n_m` adds $Cy$ to it. The
accumulator is kept as two interleaved partial sums, `r[0..n-1]` and
`t[1..n-1]`: the product $C_j y$ of limb $j$ of the cofactor goes into
limbs $j, j+1$ of `r` for even $j$ and of `t` for odd $j$, which keeps
the carry chains short and lets the multiplications for all even $j$
and all odd $j$ be issued independently. `_reduce_big_sum_n` adds the two
partial sums and subtracts the $n$-limb `limit` until the result is less
than it, leaving the reduced value in `r`.

The kernels are written as straight-line code guarded by compile-time
conditions instead of loops: the accumulators stay in registers only if
every loop indexing them is completely unrolled, which compilers do not
do reliably. For the same reason, a loop over a compile-time number of primes that calls
`_big_addmul_n_m` should be preceded by the macro `CRT_UNROLL_PRIMES`.

A Chinese remaindering for $k$ primes $p_1, \ldots, p_k$ with product
$P$ of $n$ limbs uses them as follows. With the cofactors
$C_j = P/p_j$ (of $m \le n$ limbs) and the scaled residues
$y_j = x_j (P/p_j)^{-1} \bmod p_j$, the sum $\sum_j C_j y_j$ is congruent
to the solution modulo $P$ and less than $\sum_j (P/p_j) p_j = kP$. So
`_big_mul_n_m` for the first prime, `_big_addmul_n_m` for the others and
one `_reduce_big_sum_n` give the solution in $[0, P)$, at the cost of
$km$ limb multiplications and fewer than $k$ subtractions of $P$. See
[](algorithms/integers/multimodular-reduction) for the use of such sums
at the leaves of a subproduct tree.

::: macro
CRT_UNROLL_PRIMES

A pragma asking the compiler to unroll the following loop over the
primes (`#pragma GCC unroll 8` for GCC-compatible compilers, nothing
otherwise).
:::

::: macro
CRT_DEFINE(NP, N, M)

Defines a static function with the `fft_small` CRT callback signature
that reconstructs output coefficients from `NP` prime images, where `N`
is the limb length of the product of the primes and `M` that of the
cofactors ($N = M$ or $M + 1$, $1 \le N \le 7$). Before instantiating,
the including file defines `CRT_FN(NP)` (the function name),
`CRT_Z_TYPE` (the output pointer type), `CRT_HEAD` (declarations run
once), `CRT_EMIT(zi, r, NP, N, M)` (which writes output coefficient `zi`
from the reconstructed $N$-limb value `r`; the half product of the primes
is available as `Mhalf` for signed outputs) and `CRT_TAIL` (statements
run at the end), and afterwards undefines all of them. The generated
function converts blocks of transform outputs to residues and applies the
kernels above to each coefficient. This is how `fft_small` emits integer
products, signed and unsigned polynomial coefficients and windowed
products without duplicating the reconstruction loop.
:::
