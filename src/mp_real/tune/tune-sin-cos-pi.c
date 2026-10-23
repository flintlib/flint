/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Timing of _mp_real_sin_cos_pi_ui_div_ui: the kernel on pi p/q
   against the Chebyshev iteration at each order, and against the
   complex alternative -- Newton-Taylor steps for the reciprocal q-th
   root w = exp(-i pi p/q) of (-1)^p, w0 (1 - u)^(-1/q) with
   u = 1 - (-1)^p w0^q, the powers by complex squarings and products
   by the short w0 -- written here in the same ball arithmetic.

   usage: tune-sin-cos-pi n [both [kbmax]]

   For the size n (limbs), prints the kernel's time in units of one
   multiplication M(n), then for denominators q of kb bits (q = 3 2^(kb-2)
   + 1, p = q/3 made coprime) the ratio of the Chebyshev time to the
   kernel time for the orders 4, 6, 8, 12, 16 and the tuned choice,
   and of the complex iteration (order 8).  Both outputs are computed
   unless both = 0.  The scan stops once the Chebyshev iteration is
   30% slower than the kernel. */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "flint.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "mp_real.h"

static double
wall(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

/* the minimum over three batches of at least 20 ms */
#define TIMEIT(res, stmt) do { \
    double _best = 1e100, _t0, _t1; slong _reps = 1, _k; int _r; \
    for (;;) { _t0 = wall(); for (_k = 0; _k < _reps; _k++) { stmt; } \
        _t1 = wall(); if (_t1 - _t0 > 0.02) break; _reps *= 2; } \
    for (_r = 0; _r < 3; _r++) { _t0 = wall(); for (_k = 0; _k < _reps; _k++) { stmt; } \
        _t1 = wall(); if ((_t1 - _t0) / _reps < _best) _best = (_t1 - _t0) / _reps; } \
    res = _best; } while (0)

/* the complex iteration ****************************************************/

/* (re, im) = (a + ib)^2 = (a + b)(a - b) + 2ab i */
static void
csqr(mp_real_t re, mp_real_t im, const mp_real_t a, const mp_real_t b, slong p)
{
    mp_real_t s, d;
    mp_real_init(s);
    mp_real_init(d);
    mp_real_add(s, a, b, p + 1);
    mp_real_sub(d, a, b, p + 1);
    mp_real_mul(im, a, b, p);
    mp_real_mul_2exp_si(im, im, 1);
    mp_real_mul(re, s, d, p);
    mp_real_clear(s);
    mp_real_clear(d);
}

/* one step from w0 = c0 + i s0 (exact) at p limbs: (xs, ss) = w0 (1 - u)^(-1/q) */
static void
cplx_step(mp_real_t xs, mp_real_t ss, const mp_real_t c0, const mp_real_t s0,
    ulong q, int eps, slong p)
{
    mp_real_t ar, ai, nr, ni, ur, ui, pr, pi, Wr, Wi, one, tr, ti;
    slong i, j, uu, P = FLINT_BITS * p, kb = FLINT_BIT_COUNT(q), pl, cp;

    mp_real_init(ar); mp_real_init(ai); mp_real_init(nr); mp_real_init(ni);
    mp_real_init(ur); mp_real_init(ui); mp_real_init(pr); mp_real_init(pi);
    mp_real_init(Wr); mp_real_init(Wi); mp_real_init(one); mp_real_init(tr);
    mp_real_init(ti);
    mp_real_set_ui(one, 1);
    pl = p + (kb + 8 + FLINT_BITS - 1) / FLINT_BITS;

    /* w0^q, left to right */
    mp_real_set(ar, c0);
    mp_real_set(ai, s0);
    for (i = kb - 2; i >= 0; i--)
    {
        csqr(nr, ni, ar, ai, pl);
        mp_real_swap(ar, nr);
        mp_real_swap(ai, ni);
        if ((q >> i) & 1)
            mp_real_mul_complex(ar, ai, ar, ai, c0, s0, pl);
    }

    /* u = 1 - eps w0^q */
    if (eps < 0)
    {
        mp_real_neg(ar, ar);
        mp_real_neg(ai, ai);
    }
    mp_real_sub(ur, one, ar, pl);
    mp_real_neg(ui, ai);
    uu = FLINT_MAX(mp_real_abs_bound_lt_2exp_si(ur), mp_real_abs_bound_lt_2exp_si(ui)) + 1;

    /* W = sum_{j>=1} b_j u^j, b_1 = 1/q, b_(j+1) = b_j (1 + j q) / (q (j+1)) */
    mp_real_zero(Wr);
    mp_real_zero(Wi);
    for (j = 1; j * uu - (kb - 1) >= -P - 4; j++)
    {
        cp = FLINT_MAX(2, p + (j * uu) / FLINT_BITS + 1);
        if (j == 1)
        {
            mp_real_set(pr, ur);
            mp_real_set(pi, ui);
        }
        else
            mp_real_mul_complex(pr, pi, pr, pi, ur, ui, cp);
        mp_real_set(tr, pr);
        mp_real_set(ti, pi);
        for (i = 1; i < j; i++)
        {
            mp_real_mul_ui(tr, tr, 1 + i * q, cp);
            mp_real_div_ui(tr, tr, q * (i + 1), cp);
            mp_real_mul_ui(ti, ti, 1 + i * q, cp);
            mp_real_div_ui(ti, ti, q * (i + 1), cp);
        }
        mp_real_div_ui(tr, tr, q, cp);
        mp_real_div_ui(ti, ti, q, cp);
        mp_real_add(Wr, Wr, tr, p);
        mp_real_add(Wi, Wi, ti, p);
    }

    mp_real_mul_complex(Wr, Wi, Wr, Wi, c0, s0, p);
    mp_real_add(xs, c0, Wr, p);
    mp_real_add(ss, s0, Wi, p);

    mp_real_clear(ar); mp_real_clear(ai); mp_real_clear(nr); mp_real_clear(ni);
    mp_real_clear(ur); mp_real_clear(ui); mp_real_clear(pr); mp_real_clear(pi);
    mp_real_clear(Wr); mp_real_clear(Wi); mp_real_clear(one); mp_real_clear(tr);
    mp_real_clear(ti);
}

static void
trunc_mid(mp_real_t x0, const mp_real_t xb, slong k)
{
    mp_real_set(x0, xb);
    if (x0->size > k)
    {
        flint_mpn_copyi(x0->d, x0->d + x0->size - k, k);
        x0->size = k;
    }
    x0->err = 0;
    while (x0->size > 0 && x0->d[0] == 0)
    {
        flint_mpn_copyi(x0->d, x0->d + 1, x0->size - 1);
        x0->size--;
    }
}

/* cos, sin (pi p/q) by the complex iteration of order r, from the
   kernel below 16 limbs */
static void
cplx_rec(mp_real_t xs, mp_real_t ss, ulong p, ulong q, slong n, int r)
{
    slong need, m, kb = FLINT_BIT_COUNT(q);
    mp_real_t xb, sb, x0, s0;

    need = (FLINT_BITS * n + 16 + (r - 1) * (kb + 2) + r - 1) / r;
    m = (need + FLINT_BITS - 1) / FLINT_BITS;
    if (n <= 16 || m >= n)
    {
        mp_real_sin_cos_pi_ui_div_ui(ss, xs, p, q, n);
        return;
    }

    mp_real_init(xb); mp_real_init(sb); mp_real_init(x0); mp_real_init(s0);
    cplx_rec(xb, sb, p, q, m, r);
    trunc_mid(x0, xb, m + 1);
    trunc_mid(s0, sb, m + 1);
    cplx_step(xs, ss, x0, s0, q, (p & 1) ? -1 : 1, n + 1);
    mp_real_clear(xb); mp_real_clear(sb); mp_real_clear(x0); mp_real_clear(s0);
}

/***************************************************************************/

int
main(int argc, char ** argv)
{
    slong n, kb, kbmax;
    int both, i;
    int orders[5] = {4, 6, 8, 12, 16};
    nn_ptr ys, yc;
    mp_real_t a, b, c;
    flint_rand_t state;
    double tm, tk, t;
    ulong err;

    if (argc < 2)
    {
        printf("usage: tune-sin-cos-pi n [both [kbmax]]\n");
        return 1;
    }

    n = atol(argv[1]);
    both = (argc > 2) ? atoi(argv[2]) : 1;
    kbmax = (argc > 3) ? atol(argv[3]) : FLINT_BITS - 9;

    flint_rand_init(state);
    ys = flint_malloc((n + 1) * sizeof(ulong));
    yc = flint_malloc((n + 1) * sizeof(ulong));
    mp_real_init(a);
    mp_real_init(b);
    mp_real_init(c);

    mp_real_fit_length(a, n);
    mp_real_fit_length(b, n);
    flint_mpn_urandomb(a->d, state, FLINT_BITS * n);
    flint_mpn_urandomb(b->d, state, FLINT_BITS * n);
    a->d[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
    b->d[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
    a->size = b->size = n;
    TIMEIT(tm, mp_real_mul(c, a, b, n));

    /* the kernel time does not depend on q */
    tk = 1e100;
    for (kb = 10; kb <= 50; kb += 20)
    {
        ulong q = (UWORD(3) << (kb - 2)) + 1, p = q / 3;
        while (n_gcd(p, q) != 1)
            p--;
        TIMEIT(t, _mp_real_sin_cos_pi_ui_div_ui_tune(ys, both ? yc : NULL,
            &err, p, q, n, 1, 0));
        tk = FLINT_MIN(tk, t);
    }

    printf("n = %ld, %s: M(n) = %.3g s, kernel = %.1f M(n)\n", n,
        both ? "sin and cos" : "sin only", tm, tk / tm);
    printf("  kb   Chebyshev / kernel (orders 4 6 8 12 16, tuned)   complex / kernel\n");

    for (kb = 3; kb <= kbmax; kb += (kb < 12) ? 1 : 3)
    {
        ulong q = (UWORD(3) << (kb - 2)) + 1, p = q / 3;
        double tbest = 1e100;

        while (n_gcd(p, q) != 1)
            p--;

        printf("  %2ld  ", kb);
        for (i = 0; i < 5; i++)
        {
            TIMEIT(t, _mp_real_sin_cos_pi_ui_div_ui_tune(ys, both ? yc : NULL,
                &err, p, q, n, 2, orders[i]));
            printf(" %5.2f", t / tk);
            tbest = FLINT_MIN(tbest, t);
        }
        TIMEIT(t, _mp_real_sin_cos_pi_ui_div_ui_tune(ys, both ? yc : NULL,
            &err, p, q, n, 2, 0));
        printf("  %5.2f", t / tk);

        if (both && n > 16)
        {
            TIMEIT(t, cplx_rec(a, b, p, q, n, 8));
            printf("   %5.2f", t / tk);
        }
        printf("\n");
        fflush(stdout);

        if (tbest > 1.3 * tk)
            break;
    }

    flint_free(ys);
    flint_free(yc);
    mp_real_clear(a);
    mp_real_clear(b);
    mp_real_clear(c);
    flint_rand_clear(state);
    return 0;
}
