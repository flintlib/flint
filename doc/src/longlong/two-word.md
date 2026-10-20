---
title: Two-word arithmetic with longlong.h
summary: How to form full products, accumulate sums of products in two or three words, and reduce a two-word integer modulo a word with the longlong.h macros.
order: 1
---

C has no portable integer type of twice the word size, but number
theoretic code constantly needs one: the product of two residues modulo
a 64-bit modulus has up to 128 bits. The macros of `longlong.h` give
access to the processor's full-width multiplication, add-with-carry and
double-word division instructions, with portable fallbacks.
Higher-level functions are usually preferable:
`n_mulmod2_preinv` and the [](nmod) functions for modular products, and
the [](mpn_extras) functions for arrays of limbs.

## Full products and remainders {#product}

`umul_ppmm(hi, lo, a, b)` sets `hi` and `lo` to the high and low words of
$ab$. To reduce the product modulo $d$ with `udiv_qrnnd`, the high word
must be smaller than $d$; if it may not be, reduce it first (any
$h \bmod d$ gives the same remainder). In a loop with a fixed divisor,
`udiv_qrnnd_preinv` or `n_ll_mod_preinv` replace the hardware division by
two multiplications.

## Accumulating sums of products {#accumulate}

A sum of $k$ products of words is less than $k B^2$, so for $k \le B$ it
fits in three words. Each product is added to a three-word accumulator
with `add_sssaaaaaa`, passing 0 as the top word of the product. If the
inputs are known to be smaller (for example, at most $B/2$), a two-word
accumulator may suffice for a bounded number of terms; this is the
principle of the delayed reduction in [](nmod_vec) dot products.

## Example {#example}

The program below works with 64-bit words.

```c example name=longlong_two_word
#include <flint/flint.h>
#include <flint/longlong.h>

int main(void)
{
    ulong a = UWORD(0xfedcba9876543210), b = UWORD(0x0123456789abcdef);
    ulong d = UWORD(1000000000000000003);
    ulong hi, lo, s2, s1, s0, q, r;
    ulong x[4] = {UWORD_MAX, UWORD_MAX - 1, 12345, 678};
    ulong y[4] = {UWORD_MAX, 3, 98765, 4321};
    int i;

    /* the full product a*b = hi*2^64 + lo */
    umul_ppmm(hi, lo, a, b);
    flint_printf("a*b = %wx:%016wx\n", hi, lo);

    /* (a*b) mod d; udiv_qrnnd needs its high word below d */
    udiv_qrnnd(q, r, hi % d, lo, d);
    flint_printf("a*b mod d = %wu\n", r);

    /* a dot product accumulated in three words */
    s2 = s1 = s0 = 0;
    for (i = 0; i < 4; i++)
    {
        umul_ppmm(hi, lo, x[i], y[i]);
        add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), hi, lo);
    }
    flint_printf("dot = %wx:%016wx:%016wx\n", s2, s1, s0);

    /* signed product and division, rounding towards zero */
    smul_ppmm(hi, lo, WORD(-3), WORD(5));
    sdiv_qrnnd(q, r, hi, lo, WORD(4));
    flint_printf("-15 / 4: q = %wd, r = %wd\n", (slong) q, (slong) r);

    flint_printf("bits(1000) = %wu, clz(1) = %d, ctz(40) = %d\n",
                 FLINT_BIT_COUNT(1000), flint_clz(UWORD(1)), flint_ctz(UWORD(40)));
    flint_printf("flog2(1000) = %wu, clog2(1000) = %wu, clog2(1024) = %wu\n",
                 FLINT_FLOG2(1000), FLINT_CLOG2(1000), FLINT_CLOG2(1024));
    return 0;
}
```

```output
a*b = 121fa00ad77d742:2236d88fe5618cf0
a*b mod d = 105526465657053906
dot = 1:0000000000000001:0000000048d90a86
-15 / 4: q = -3, r = -3
bits(1000) = 10, clz(1) = 63, ctz(40) = 3
flog2(1000) = 9, clog2(1000) = 10, clog2(1024) = 10
```
