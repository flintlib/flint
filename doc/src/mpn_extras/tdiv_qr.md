---
title: Euclidean division of limb arrays
functions: flint_mpn_tdiv_qr flint_mpn_cdiv_qr flint_mpn_ndiv_qr flint_mpn_divapprox flint_mpn_inv
---

`flint_mpn_tdiv_qr` replaces GMP's `mpn_tdiv_qr`, with the same interface
and aliasing rules, and the other functions below are built on the
same machinery. The division is done entirely by FLINT (except for GMP's
assembly `mpn_divrem_1` and, where it helps, `mpn_divrem_2`), which makes
it faster than GMP's for all shapes, and asymptotically faster for huge
operands, where it uses Newton iteration on top of the `fft_small`
multiplication.

::: summary algorithms/integers/division :::

## Algorithm selection {#selection}

With $m = \mathit{an} - \mathit{bn} + 1$ quotient limbs, `_flint_mpn_tdiv_qr`
chooses, in this order:

1. **Register-based division** (`_flint_mpn_tdiv_qr_small`) for
   divisors of at most 3 limbs, and for divisors of 4 to 7 limbs with
   short dividends ($\mathit{an} \le 3\,\mathit{bn}$ when the remainder
   is wanted, $\mathit{an} \le \mathit{bn} + 2$ otherwise). The divisor and
   the partial remainder are kept in registers and each quotient limb
   takes one 3/2 division step [@GraMol2010].
2. **Newton division** when the divisor and the quotient both have at
   least `FLINT_MPN_TDIV_QR_NEWTON_CUTOFF` limbs (950 on x86-64), or the
   divisor has at least `FLINT_MPN_TDIV_QR_NEWTON_LONG_CUTOFF` limbs (608)
   and the quotient is at least twice as long: `_flint_mpn_tdiv_qr_newton`
   computes the quotient as a fixed-point number by Karp–Markstein
   iteration [@KM1997] (`_mp_real_div_newton`) and certifies its integer
   part, or `_flint_mpn_tdiv_qr_unbalanced` if $\mathit{an} \ge
   3\,\mathit{bn}$.
3. **Blockwise division with a shared inverse**
   (`_flint_mpn_tdiv_qr_unbalanced`) for $\mathit{an} \ge 4\,\mathit{bn}$
   with $\mathit{bn} \ge 64$, or $\mathit{an} \ge 3\,\mathit{bn}$ with
   $\mathit{bn} \ge 512$: the dividend is processed in blocks of
   $2\,\mathit{bn} \times \mathit{bn}$ divisions that all use one
   approximate inverse of $b$.
4. **Division with an $n$-limb precomputed inverse**
   (`_flint_mpn_tdiv_qr_preinvn`) for short divisors and long dividends:
   $\mathit{bn} \ge 32$ and $\mathit{an} \ge 4\,\mathit{bn}$, or
   $\mathit{bn} \ge 4$ and $\mathit{an} \ge 32\,\mathit{bn}$.
5. Otherwise **schoolbook and divide and conquer division**
   (`_flint_mpn_tdiv_qr_divconquer`, or `_flint_mpn_tdiv_q_divconquer`
   without the remainder): ports of GMP's code [@BurZie1998] with FLINT's
   multiplication, and for quotients about as short as the divisor, a
   short division in the style of Mulders that only forms high products
   [@HZ2011].

All cutoffs are machine dependent (`flint-mparam.h`, measured by
`src/mpn_extras/tune/tune-div.c`); the values above are those for
x86-64.

Where hardware division is fast (`FLINT_PREINVERT_LIMB_USE_NATIVE`, for
example on AMD Zen 3 and later), limb inverses are computed with it, and
short dividends use chains of hardware divisions instead of inverses:
one-limb divisors below `FLINT_MPN_DIVREM_1_HW_CUTOFF` dividend limbs
(`FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF` for normalised divisors), two-limb
divisors below `FLINT_MPN_DIV_2_HW_CUTOFF` limbs (a 3/2 step made of one
hardware 2/1 division and a correction by the low divisor limb, as in
Knuth's algorithm D, with normalisation on the fly), and 3- to 7-limb
divisors for quotients shorter than `FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF`.
Two-limb divisors use GMP's assembly `mpn_divrem_2` (on a shifted copy of
the dividend) from `FLINT_MPN_DIV_2_GMP_CUTOFF` dividend limbs where it
is faster.

## Certifying a Newton quotient {#newton}

Newton division computes the quotient $a/b$ as a fixed-point number with
$m$ integer limbs and two fraction limbs, with an error below $4B^{-2}$ in
absolute value. If the first fraction limb $f$ satisfies $2 \le f \le
B - 2$, the true quotient $a/b$ lies strictly between the computed integer
part $q$ and $q + 1$, so $q = \lfloor a/b \rfloor$; then the remainder
$a - qb$ is less than $b$ and can be computed from a low product of
$\mathit{bn}$ limbs. Otherwise (with probability about $4/B$ for random
inputs) the candidate is corrected by $O(1)$ steps using a full product.

## The variants {#variants}

- `flint_mpn_tdiv_q` and `flint_mpn_tdiv_r` compute only the quotient or
  only the remainder (the latter with a temporary quotient). Without
  the remainder, several of the algorithms above save work: the short
  division forms only high products, and the Newton division skips the
  remainder unless the certification fails.
- `flint_mpn_cdiv_qr` and `flint_mpn_ndiv_qr` round the quotient up or to
  nearest (ties to even), with one extra quotient limb.
- `flint_mpn_divapprox` returns $\lfloor a/b \rfloor$ or one more,
  skipping all corrections; `flint_mpn_divapprox_fraction` adds fraction
  limbs. These are what floating-point division needs.
- `flint_mpn_inv` and `flint_mpn_invapprox` compute
  $\lfloor B^n / x \rfloor$ (or one more) without forming the numerator.
- `flint_mpn_div` and `flint_mpn_divexact` are the exact divisions; the
  latter uses Hensel division (see [](mpn_extras#hensel)).

## Example {#example}

```c example name=mpn_tdiv_qr_basic
#include <gmp.h>
#include <flint/flint.h>
#include <flint/mpn_extras.h>

int main(void)
{
    flint_rand_t state;
    mp_size_t an = 100, bn = 37, qn = an - bn + 1;
    mp_limb_t a[100], b[37], q[64], r[37], q2[64], r2[37], t[101];

    flint_rand_init(state);
    flint_mpn_urandomb(a, state, an * FLINT_BITS);
    flint_mpn_urandomb(b, state, bn * FLINT_BITS);
    b[bn - 1] |= 1;     /* nonzero top limb */

    flint_mpn_tdiv_qr(q, r, a, an, b, bn);

    /* compare with GMP */
    mpn_tdiv_qr(q2, r2, 0, a, an, b, bn);
    flint_printf("agrees with GMP: %d\n",
        mpn_cmp(q, q2, qn) == 0 && mpn_cmp(r, r2, bn) == 0);

    /* check a = q b + r and r < b */
    flint_mpn_mul(t, q, qn, b, bn);
    mpn_add(t, t, an + 1, r, bn);
    flint_printf("a = q b + r: %d, r < b: %d\n",
        mpn_cmp(t, a, an) == 0 && t[an] == 0, mpn_cmp(r, b, bn) < 0);

    /* approximate quotient: q or q + 1 */
    flint_mpn_divapprox(q2, a, an, b, bn);
    mpn_sub_n(q2, q2, q, qn);
    flint_printf("divapprox - tdiv_q <= 1: %d\n",
        q2[0] <= 1 && flint_mpn_zero_p(q2 + 1, qn - 1));

    flint_rand_clear(state);
    return 0;
}
```
```output
agrees with GMP: 1
a = q b + r: 1, r < b: 1
divapprox - tdiv_q <= 1: 1
```
