/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Timing of arf_*_approx (fixed_*_newton) against the previous arb Newton
   iterations in arf arithmetic and against the plain arf (MPFR/GMP)
   functions. */

#include "arb.h"
#include "profiler.h"

#define GUARD_BITS 32
#define INV_NEWTON_CUTOFF 24000
#define RSQRT_NEWTON_CUTOFF 4000

static void
_old_inv_newton(arf_t res, const arf_t x, slong prec)
{
    slong wp = prec + GUARD_BITS;
    slong hp = prec / 2 + GUARD_BITS;

    if (prec < INV_NEWTON_CUTOFF)
    {
        arf_set_round(res, x, wp, ARF_RND_DOWN);
        arf_ui_div(res, 1, res, wp, ARF_RND_DOWN);
    }
    else
    {
        arf_t r, t;

        arf_init(r);
        arf_init(t);

        _old_inv_newton(r, x, hp);

        /* r - r*(x*r - 1) */

        if (arf_bits(x) <= wp)
        {
            arf_mul(t, x, r, wp, ARF_RND_DOWN);
        }
        else
        {
            arf_set_round(t, x, wp, ARF_RND_DOWN);
            arf_mul(t, t, r, wp, ARF_RND_DOWN);
        }

        arf_sub_ui(t, t, 1, hp, ARF_RND_DOWN);
        arf_mul(t, t, r, hp, ARF_RND_DOWN);
        arf_sub(res, r, t, wp, ARF_RND_DOWN);

        arf_clear(r);
        arf_clear(t);
    }
}

/* Karp-Markstein */
static void
_old_div_newton(arf_t res, const arf_t x, const arf_t y, slong prec)
{
    arf_t xn, yn, t;

    slong wp = prec + GUARD_BITS;
    slong hp = prec / 2 + GUARD_BITS;

    arf_init(xn);
    arf_init(yn);
    arf_init(t);

    _old_inv_newton(xn, y, hp);
    arf_set_round(t, x, hp, ARF_RND_DOWN);
    arf_mul(yn, xn, t, hp, ARF_RND_DOWN);
    arf_mul(t, y, yn, wp, ARF_RND_DOWN);
    arf_sub(t, x, t, hp, ARF_RND_DOWN);
    arf_mul(t, t, xn, hp, ARF_RND_DOWN);
    arf_add(res, yn, t, wp, ARF_RND_DOWN);

    arf_clear(xn);
    arf_clear(yn);
    arf_clear(t);
}
static void
_old_rsqrt_newton(arf_t res, const arf_t x, slong prec)
{
    slong wp = prec + GUARD_BITS;
    slong hp = prec / 2 + GUARD_BITS;

    if (prec < RSQRT_NEWTON_CUTOFF)
    {
        arf_set_round(res, x, wp, ARF_RND_DOWN);
        arf_rsqrt(res, res, wp, ARF_RND_DOWN);
    }
    else
    {
        arf_t r, t, u;

        arf_init(r);
        arf_init(t);
        arf_init(u);

        _old_rsqrt_newton(r, x, hp);

        /* r - r*(x*r^2 - 1)/2 */

        arf_mul(t, r, r, wp, ARF_RND_DOWN);

        if (arf_bits(x) <= wp)
        {
            arf_mul(t, t, x, wp, ARF_RND_DOWN);
        }
        else
        {
            arf_set_round(u, x, wp, ARF_RND_DOWN);
            arf_mul(t, t, u, wp, ARF_RND_DOWN);
        }

        arf_sub_ui(t, t, 1, hp, ARF_RND_DOWN);
        arf_mul_2exp_si(t, t, -1);
        arf_mul(t, t, r, hp, ARF_RND_DOWN);

        arf_sub(res, r, t, wp, ARF_RND_DOWN);

        arf_clear(r);
        arf_clear(t);
        arf_clear(u);
    }
}

static void
_old_sqrt_newton(arf_t res, const arf_t x, slong prec)
{
    arf_t t, u, v;

    slong wp = prec + GUARD_BITS;
    slong hp = prec / 2 + GUARD_BITS;

    arf_init(t);
    arf_init(u);
    arf_init(v);

    _old_rsqrt_newton(t, x, hp);

    if (arf_bits(x) <= hp)
    {
        arf_mul(v, t, x, hp, ARF_RND_DOWN);
    }
    else
    {
        arf_set_round(u, x, hp, ARF_RND_DOWN);
        arf_mul(v, t, u, hp, ARF_RND_DOWN);
    }

    arf_mul(u, v, v, wp, ARF_RND_DOWN);
    arf_sub(u, x, u, hp, ARF_RND_DOWN);
    arf_mul(u, u, t, wp, ARF_RND_DOWN);
    arf_mul_2exp_si(u, u, -1);
    arf_add(res, v, u, wp, ARF_RND_DOWN);

    arf_clear(t);
    arf_clear(u);
    arf_clear(v);
}

#define TIME(expr, res) \
    do { timeit_t __t; slong __r; expr; \
         TIMEIT_REPEAT(__t, __r) { expr; } TIMEIT_END_REPEAT(__t, __r); \
         res = (double) __t->wall * 0.001 / __r; } while (0)

int main(void)
{
    flint_rand_t state;
    slong prec;
    arf_t x, y, r;
    double t1, t2, t3;

    flint_rand_init(state);
    arf_init(x); arf_init(y); arf_init(r);

    flint_printf("%10s | %10s %10s %10s | %10s %10s %10s\n", "prec", "arf_div", "old_div", "div_fast", "arf_inv", "old_inv", "inv_fast");
    for (prec = 4000; prec <= 8000000; prec *= 2)
    {
        { fmpz_t f; fmpz_init(f); fmpz_randbits(f, state, prec); fmpz_abs(f, f); fmpz_setbit(f, prec - 1); arf_set_fmpz(x, f);
          fmpz_randbits(f, state, prec); fmpz_abs(f, f); fmpz_setbit(f, prec - 1); arf_set_fmpz(y, f); fmpz_clear(f); }
        TIME(arf_div(r, x, y, prec, ARF_RND_DOWN), t1);
        TIME(_old_div_newton(r, x, y, prec), t2);
        TIME(arf_div(r, x, y, prec, ARF_RND_FAST), t3);
        flint_printf("%10wd | %10.3e %10.3e %10.3e |", prec, t1, t2, t3);
        TIME(arf_ui_div(r, 1, y, prec, ARF_RND_DOWN), t1);
        TIME(_old_inv_newton(r, y, prec), t2);
        TIME(arf_ui_div(r, 1, y, prec, ARF_RND_FAST), t3);
        flint_printf(" %10.3e %10.3e %10.3e\n", t1, t2, t3);
    }

    flint_printf("%10s | %10s %10s %10s | %10s %10s %10s\n", "prec", "arf_sqrt", "old_sqrt", "sqrt_fast", "arf_rsqrt", "old_rsqrt", "rsqrt_fast");
    for (prec = 4000; prec <= 8000000; prec *= 2)
    {
        { fmpz_t f; fmpz_init(f); fmpz_randbits(f, state, prec); fmpz_abs(f, f); fmpz_setbit(f, prec - 1); arf_set_fmpz(x, f); fmpz_clear(f); }
        TIME(arf_sqrt(r, x, prec, ARF_RND_DOWN), t1);
        TIME(_old_sqrt_newton(r, x, prec), t2);
        TIME(arf_sqrt(r, x, prec, ARF_RND_FAST), t3);
        flint_printf("%10wd | %10.3e %10.3e %10.3e |", prec, t1, t2, t3);
        TIME(arf_rsqrt(r, x, prec, ARF_RND_DOWN), t1);
        TIME(_old_rsqrt_newton(r, x, prec), t2);
        TIME(arf_rsqrt(r, x, prec, ARF_RND_FAST), t3);
        flint_printf(" %10.3e %10.3e %10.3e\n", t1, t2, t3);
    }

    arf_clear(x); arf_clear(y); arf_clear(r);
    flint_rand_clear(state);
    return 0;
}
