---
title: Mathematical constants
functions: mp_real_const_pi4 _mp_real_const_pi4 _mp_real_const_clear_cache mp_real_const_func _mp_real_const_arb _mp_real_const_euler_tune
---

`mp_real` computes twelve constants: $\pi/4$, $\log 2$, Euler's constant
$\gamma$, $e$, $\log 10$, Catalan's constant $G$, $\zeta(3)$,
$\zeta(5)$, $\Gamma(1/3)$, $\Gamma(1/4)$, $2/\pi$ and $1/\log 2 - 1$.
They are the implementations behind `arb_const_pi`, `arb_const_log2`,
`arb_const_euler`, `arb_const_e`, `arb_const_log10`,
`arb_const_catalan`, `arb_const_apery`, `arb_const_zeta5` and the
$\Gamma$ values of `arb_hypgeom_gamma_fmpq`, and they are used inside
the module by the argument reductions of the elementary functions.

## Two forms {#forms}

The ball functions `mp_real_const_*` return a ball for the constant
accurate to about $n$ limbs. The limb functions `_mp_real_const_*`
return the exact floor $\lfloor c B^n \rfloor$ as a fixed-point number
with $n$ fraction limbs, plus a units limb for the constants above 1,
and set `err` to 1: the true value lies in $[y, y + 1)$ ulps. An exact
floor is a canonical value, independent of how it was computed. It is
determined from a ball computed with a few guard limbs, recomputed with
more guard limbs in the rare case that the ball contains a multiple of
$B^{-n}$.

## Caching {#cache}

With `cache` = 0 the constant is computed from scratch and nothing is
kept. With `cache` = 1 the value comes from a process-wide cache that
holds the floor at the largest precision requested so far. Floors nest:
the top limbs of $\lfloor cB^N \rfloor$ are $\lfloor cB^n \rfloor$ for
$n \le N$, so a shorter request is a copy of the top limbs. A longer
request recomputes the entry at $\max(n + 5, 1.5N)$ limbs, so that a
sequence of growing requests costs a bounded factor more than the
largest. All threads read the same entry; one that needs a longer entry
builds it under the cache's lock, and readers never wait. The ball
version with `cache` = 1 returns the cached floor with a radius of one
ulp. $\pi/4$, $\log 2$, $2/\pi$ and $1/\log 2 - 1$ also come from static
tables of 64 limbs, so that the elementary functions at small precision
touch no cache at all. `_mp_real_const_clear_cache` (or
`flint_cleanup_master`) frees the caches.

The `arb` constants call `_mp_real_const_arb`, which computes with
`cache` = 0 and leaves caching to `arb`'s own caches.

## Algorithms {#algorithms}

All constants except $\log 10$ are hypergeometric series summed by
`mp_real_hypgeom_series`, with formulas chosen from y-cruncher's
formula collection as the fastest in this implementation:

| Constant | Formula |
|---|---|
| $\pi$ | Chudnovsky's series [@CC1988], 47 bits per term |
| $\log 2$ | Zuniga's series [@Zun2025], 11.9 bits per term |
| $e$ | $\sum 1/k!$ |
| $\gamma$ | Brent–McMillan [@BM1980] with the bound of [@BJ2013], two splittings in dual numbers |
| $\log 10$ | $\log 2 + \log 5$ from a three-term Machin-type formula |
| $G$ | Pilehrood's series [@PP2010] |
| $\zeta(3)$ | a series of Zuniga (2023), 2.05 bits per term |
| $\zeta(5)$ | Sun's identity $\zeta(5) = (3S + 56\pi^2\zeta(3))/540$ |
| $\Gamma(1/3)$ | $(810^{1/4}\pi/X)^{1/3}$ with Guillera's series $X$ |
| $\Gamma(1/4)$ | Ebisu's series, or $\sqrt{(2\pi)^{3/2}/\operatorname{agm}(1, \sqrt 2)}$ from $2^{16}$ limbs |

The independent parts (the series and $\pi$ in $\zeta(5)$ and the
$\Gamma$ values, the three logarithms of $\log 10$, the logarithms
needed by $\gamma$) run on separate threads. The details, including the
dual-number splitting for $\gamma$, are in
[](algorithms/floating-point/constants).

::: summary algorithms/floating-point/constants :::

## Example {#example}

```c example name=mp_real_constants
#include <flint/flint.h>
#include <flint/mp_real.h>
#include <flint/arb.h>

int main(void)
{
    mp_real_t c;
    arb_t a;
    ulong f[3], err;
    slong n = mp_real_prec_bits(200);

    mp_real_init(c);
    arb_init(a);

    /* zeta(3) as a ball, computed from scratch */
    mp_real_const_zeta3(c, n, 0);
    mp_real_get_arb(a, c);
    flint_printf("zeta(3) = ");
    arb_printn(a, 50, 0);
    flint_printf("\n");

    /* the cached floor of e B^2: two fraction limbs and a units limb */
    _mp_real_const_e(f, &err, 2, 1);
    flint_printf("floor(e B^2) = %wu %016wx %016wx, err = %wu\n",
        f[2], f[1], f[0], err);

    /* Catalan's constant through the arb entry point */
    _mp_real_const_arb(a, mp_real_const_catalan, 128);
    flint_printf("G = ");
    arb_printn(a, 30, 0);
    flint_printf("\n");

    mp_real_clear(c);
    arb_clear(a);
    flint_cleanup_master();
    return 0;
}
```

```output
zeta(3) = [1.2020569031595942853997381615114499907649862923405 +/- 1.12e-51]
floor(e B^2) = 2 b7e151628aed2a6a bf7158809cf4f3c7, err = 1
G = [0.915965594177219015054603514932 +/- 3.85e-31]
```

The example assumes 64-bit limbs; the hexadecimal digits of $e$ are
$2.\mathrm{B7E151628AED2A6ABF7158809CF4F3C7}\ldots$
