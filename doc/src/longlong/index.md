---
title: longlong.h
header: longlong.h
description: Macros for double-word and multi-word arithmetic on limbs, bit counting and division with precomputed inverses.
category: Foundations
order: 6
---

This header provides the primitive operations on machine words that
multiprecision arithmetic is built from: counting leading and trailing
zero bits, adding and subtracting integers of two to eleven limbs held in
separate variables, the full two-word product of two words, and the
division of a two-word integer by a word. Most of them are macros that
assign to their first arguments, in the style of GMP's `longlong.h`, from
which the names come: in `add_ssaaaa(s1, s0, a1, a0, b1, b0)` each letter
stands for one argument, `s` for a word of the sum and `a`, `b` for the
words of the operands, most significant first.

The macros are implemented with inline assembly on x86-64 and ARM64 (and
with intrinsics under MSVC), and in portable C elsewhere. Throughout,
$B = 2^{\bits}$, and all words are `ulong` values. The guide
[](longlong/two-word) shows how they are combined.

The arguments of the macros may be evaluated more than once, so they
should be variables or expressions without side effects. With the
portable C versions, an output may be the same variable as the input
word of the same position (`s0` as `a0` or `b0`, and so on), but not as an
input word of a different position. For speed, write sums in place
(`add_ssaaaa(s1, s0, s1, s0, b1, b0)`) and put constants in the `b`
arguments of additions and the `s` arguments of subtractions.

## Bit counting {#bits}

::: macro
flint_clz(x)
flint_ctz(x)

The number of leading (most significant) and trailing (least
significant) zero bits of the word $x$, as an `int`. The result is
undefined for $x = 0$. `flint_clz(x)` is the shift that normalises $x$,
that is, makes its top bit 1.
:::

::: function
flint_bitcnt_t FLINT_BIT_COUNT(ulong x)

Returns the number of bits of $x$, that is, $\lfloor \log_2 x \rfloor + 1$
for $x \ge 1$, and 0 for $x = 0$. The result is unsigned, so expressions
such as `FLINT_BIT_COUNT(x) - 1` wrap around for $x = 0$.
:::

::: macro
FLINT_FLOG2(x)
FLINT_CLOG2(x)

$\lfloor \log_2 x \rfloor$ and $\lceil \log_2 x \rceil$, for $x \ge 1$.
:::

## Addition and subtraction {#addsub}

All these macros work modulo a power of $B$: the carry out of the most
significant word, or the borrow, is discarded.

::: macro
add_ssaaaa(s1, s0, a1, a0, b1, b0)
sub_ddmmss(d1, d0, m1, m0, s1, s0)

Set $(s_1, s_0)$ to $(a_1 B + a_0) + (b_1 B + b_0) \bmod B^2$, and
$(d_1, d_0)$ to $(m_1 B + m_0) - (s_1 B + s_0) \bmod B^2$.
:::

::: macro
add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0)
sub_dddmmmsss(d2, d1, d0, m2, m1, m0, s2, s1, s0)

The same for three-word integers, modulo $B^3$.
:::

::: macro
add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)
sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)
add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0)
sub_dddddmmmmmsssss(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0)
add_ssssssaaaaaaaaaaaa(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0)
sub_ddddddmmmmmmssssss(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0)
add_sssssssaaaaaaaaaaaaaa(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0)
sub_dddddddmmmmmmmsssssss(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0)
add_ssssssssaaaaaaaaaaaaaaaa(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0)
sub_ddddddddmmmmmmmmssssssss(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0)
add_sssssssssaaaaaaaaaaaaaaaaaa(s8, s7, s6, s5, s4, s3, s2, s1, s0, a8, a7, a6, a5, a4, a3, a2, a1, a0, b8, b7, b6, b5, b4, b3, b2, b1, b0)
sub_dddddddddmmmmmmmmmsssssssss(s8, s7, s6, s5, s4, s3, s2, s1, s0, a8, a7, a6, a5, a4, a3, a2, a1, a0, b8, b7, b6, b5, b4, b3, b2, b1, b0)
add_ssssssssssaaaaaaaaaaaaaaaaaaaa(s9, s8, s7, s6, s5, s4, s3, s2, s1, s0, a9, a8, a7, a6, a5, a4, a3, a2, a1, a0, b9, b8, b7, b6, b5, b4, b3, b2, b1, b0)
sub_ddddddddddmmmmmmmmmmssssssssss(s9, s8, s7, s6, s5, s4, s3, s2, s1, s0, a9, a8, a7, a6, a5, a4, a3, a2, a1, a0, b9, b8, b7, b6, b5, b4, b3, b2, b1, b0)
add_sssssssssssaaaaaaaaaaaaaaaaaaaaaa(s10, s9, s8, s7, s6, s5, s4, s3, s2, s1, s0, a10, a9, a8, a7, a6, a5, a4, a3, a2, a1, a0, b10, b9, b8, b7, b6, b5, b4, b3, b2, b1, b0)
sub_dddddddddddmmmmmmmmmmmsssssssssss(s10, s9, s8, s7, s6, s5, s4, s3, s2, s1, s0, a10, a9, a8, a7, a6, a5, a4, a3, a2, a1, a0, b10, b9, b8, b7, b6, b5, b4, b3, b2, b1, b0)

Addition and subtraction of integers of 4 to 11 words, modulo $B^4$ to
$B^{11}$.
:::

::: macro
MPN_INCR_U(ptr, size, incr)
MPN_DECR_U(ptr, size, incr)

Add or subtract the word `incr` to or from the integer of `size` limbs
at `ptr`, in place, using GMP's `mpn_add_1` and `mpn_sub_1` (so `gmp.h`
must be included). The carry or borrow out must be zero; this is checked
when assertions are enabled.
:::

## Multiplication {#multiplication}

::: macro
umul_ppmm(p1, p0, u, v)

Sets $(p_1, p_0)$ to the two-word product $p_1 B + p_0 = u v$ of the
words $u$ and $v$.
:::

::: macro
smul_ppmm(p1, p0, u, v)

Sets $(p_1, p_0)$ to the two-word product of $u$ and $v$ interpreted as
signed words, in two's complement: $p_1$ read as an `slong` is the
signed high word.
:::

## Division {#division}

::: macro
udiv_qrnnd(q, r, n1, n0, d)

Sets $q$ and $r$ to the quotient and remainder of $n_1 B + n_0$ divided
by $d$, so that $n_1 B + n_0 = q d + r$ with $0 \le r < d$. Requires
$n_1 < d$, which ensures that $q < B$; otherwise the result is
undefined, and on x86-64 the hardware division traps.
:::

::: macro
sdiv_qrnnd(q, r, n1, n0, d)

The signed version of `udiv_qrnnd`: $n_1 B + n_0$ and $d$ are signed
two's complement integers, the quotient is rounded towards zero and the
remainder has the sign of the dividend (as with the C operators `/` and
`%`). Requires that the quotient fits in an `slong`.
:::

::: macro
udiv_qrnnd_preinv(q, r, n1, n0, d, dinv)

As `udiv_qrnnd`, but faster: uses the precomputed inverse
`dinv = n_preinvert_limb_prenorm(d)` instead of a hardware division. The
divisor $d$ must be normalised (its top bit set) and $n_1 < d$. This is
the division algorithm of Möller and Granlund
[](algorithms/arithmetic/moller-granlund).
:::

## Internal macros {#internal}

::: macro
udiv_qrnnd_int(q, r, n1, n0, d)

The portable C implementation of `udiv_qrnnd` for a normalised divisor,
by schoolbook division with half words. Only defined when no assembly
version of `udiv_qrnnd` is available.
:::

::: macro
__ll_B
__ll_lowpart(t)
__ll_highpart(t)

$2^{\bits/2}$ and the low and high halves of a word, used by the portable
versions of `umul_ppmm` and `udiv_qrnnd`.
:::
