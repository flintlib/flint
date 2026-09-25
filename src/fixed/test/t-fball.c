/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "arb.h"
#include "fixed.h"

/* the fball ball arithmetic: every operation's output ball must
   contain the exact result of a random true point drawn inside each
   operand ball (with aliasing variants), the error normalization
   invariant must hold afterward, and the two series constants must
   overlap arb's at a radius within a couple of guard limbs */

/* random normalized fball; if allow_err, a random radius (an integer
   count of ulps of the bottom limb, which zero padding may place well
   below the significant limbs) */
static void
fball_randtest(fball_t x, flint_rand_t state, slong maxsize, int allow_err)
{
    slong size = n_randint(state, maxsize + 1);

    if (size == 0 && n_randint(state, 2))
    {
        fball_zero(x);
        x->exp = (slong) n_randint(state, 13) - 6;
        if (allow_err && n_randint(state, 2))
            x->err = 1 + n_randint(state, 1000);
        return;
    }

    size = FLINT_MAX(size, 1);
    fball_fit(x, size);
    flint_mpn_rrandom(x->d, state, size);
    x->d[size - 1] |= (UWORD(1) << (FLINT_BITS - 1 -
        n_randint(state, FLINT_BITS - 1)));
    x->size = size;
    x->negative = (int) n_randint(state, 2);
    x->exp = size + (slong) n_randint(state, 13) - 6;
    x->err = 0;
    if (allow_err && n_randint(state, 2))
    {
        x->err = 1 + n_randint(state, 1000);
        if (n_randint(state, 2))
            x->err = n_randtest(state) | 1;
        if (n_randint(state, 3) == 0)
        {
            /* zero padding: a radius far below the significant limbs */
            slong pad = 1 + n_randint(state, 4);
            fball_fit(x, x->size + pad);
            memmove(x->d + pad, x->d, x->size * sizeof(ulong));
            flint_mpn_zero(x->d, pad);
            x->size += pad;
        }
    }
    if (x->err == 0 && n_randint(state, 2))
    {
        /* exercise the low-zero-limb stripping */
        x->d[0] &= ~(ulong) n_randint(state, 2);
    }
    /* keep the normalization invariants */
    if (x->err == 0)
    {
        slong t = 0;
        while (t < x->size && x->d[t] == 0)
            t++;
        if (t == x->size)
        {
            fball_zero(x);
            return;
        }
        if (t > 0)
        {
            flint_mpn_copyi(x->d, x->d + t, x->size - t);
            x->size -= t;
        }
    }
}

/* a random true value inside the ball, as an exact arf */
static void
fball_random_point(arf_t t, const fball_t x, flint_rand_t state)
{
    arb_t b;
    arf_t u;

    arb_init(b);
    arf_init(u);
    fball_get_arb(b, x);

    arf_set(t, arb_midref(b));
    if (x->err != 0)
    {
        /* mid + (r / 2^30) * err * ulp, r in [-2^30, 2^30] */
        slong r = (slong) n_randint(state, UWORD(1) << 31)
                    - (slong) (UWORD(1) << 30);
        arf_set_ui(u, x->err);
        arf_mul_si(u, u, r, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul_2exp_si(u, u,
            FLINT_BITS * ((x->size == 0) ? x->exp : x->exp - x->size) - 30);
        arf_add(t, t, u, ARF_PREC_EXACT, ARF_RND_NEAR);
    }

    arb_clear(b);
    arf_clear(u);
}

static void
check_contains(const fball_t res, const arf_t truth, const char * op,
    slong iter)
{
    arb_t r;
    arb_init(r);
    fball_get_arb(r, res);
    if (!arb_contains_arf(r, truth))
    {
        flint_printf("FAIL: %s (iter %wd)\nres = ", op, iter);
        fball_print(res);
        flint_printf("truth = "); arf_printd(truth, 30);
        flint_printf("\n");
        flint_abort();
    }
    arb_clear(r);
}

static void
test_ops(flint_rand_t state, slong iters, slong maxsize, slong maxn)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        fball_t a, b, r;
        arf_t ta, tb, tr;
        slong n = 1 + n_randint(state, maxn);
        int aliased = (int) n_randint(state, 3);

        fball_init(a); fball_init(b); fball_init(r);
        arf_init(ta); arf_init(tb); arf_init(tr);

        fball_randtest(a, state, maxsize, 1);
        fball_randtest(b, state, maxsize, 1);
        fball_random_point(ta, a, state);
        fball_random_point(tb, b, state);

        /* add */
        if (aliased == 1) { fball_set(r, a); fball_add(r, r, b, n); }
        else if (aliased == 2) { fball_set(r, b); fball_add(r, a, r, n); }
        else fball_add(r, a, b, n);
        arf_add(tr, ta, tb, ARF_PREC_EXACT, ARF_RND_NEAR);
        check_contains(r, tr, "add", iter);

        /* sub (aliased variants exercise the in-place window
           build and negation) */
        if (aliased == 1) { fball_set(r, a); fball_sub(r, r, b, n); }
        else if (aliased == 2) { fball_set(r, b); fball_sub(r, a, r, n); }
        else fball_sub(r, a, b, n);
        arf_sub(tr, ta, tb, ARF_PREC_EXACT, ARF_RND_NEAR);
        check_contains(r, tr, "sub", iter);

        /* addmul_ui / submul_ui, all aliasings; b is often placed
           inside a's window (the in-place path) */
        {
            ulong c = n_randint(state, 4) ? n_randtest(state) : n_randint(state, 10);
            fball_t bb;
            arf_t tbb, u;
            int sub = n_randint(state, 2);
            fball_init(bb); arf_init(tbb); arf_init(u);
            fball_set(bb, b);
            if (n_randint(state, 2) && a->size > 2 && b->size > 0)
            {
                /* shift b to sit within a's limbs, or to reach below
                   a's bottom while ending under its top */
                fball_mul_2exp_si(bb, FLINT_BITS * ((a->exp - 2 - (slong) n_randint(state, 2)) - bb->exp));
                if (n_randint(state, 2))
                    fball_mul_2exp_si(bb, -FLINT_BITS * (slong) n_randint(state, 3));
            }
            fball_random_point(tbb, bb, state);
            arf_mul_ui(u, tbb, c, ARF_PREC_EXACT, ARF_RND_NEAR);
            if (sub)
                arf_sub(tr, ta, u, ARF_PREC_EXACT, ARF_RND_NEAR);
            else
                arf_add(tr, ta, u, ARF_PREC_EXACT, ARF_RND_NEAR);
            if (aliased == 1) { fball_set(r, a); if (sub) fball_submul_ui(r, r, bb, c, n); else fball_addmul_ui(r, r, bb, c, n); }
            else if (aliased == 2) { fball_set(r, bb); if (sub) fball_submul_ui(r, a, r, c, n); else fball_addmul_ui(r, a, r, c, n); }
            else { if (sub) fball_submul_ui(r, a, bb, c, n); else fball_addmul_ui(r, a, bb, c, n); }
            check_contains(r, tr, sub ? "submul_ui" : "addmul_ui", iter);
            /* exact inputs fitting the precision: exact result */
            if (a->err == 0 && bb->err == 0 && r->err != 0 && a->size + 2 < n
                && bb->exp - bb->size >= a->exp - a->size && bb->exp + 1 <= a->exp)
            {
                flint_printf("FAIL: addmul_ui inexact (iter %wd)\n", iter);
                fball_print(a); fball_print(bb); fball_print(r);
                flint_abort();
            }
            fball_clear(bb); arf_clear(tbb); arf_clear(u);
        }

        /* the same ball twice (doubling and cancellation), also in
           place */
        if (aliased == 1) { fball_set(r, a); fball_add(r, r, r, n); }
        else fball_add(r, a, a, n);
        arf_add(tr, ta, ta, ARF_PREC_EXACT, ARF_RND_NEAR);
        check_contains(r, tr, "add a a", iter);
        if (aliased == 1) { fball_set(r, a); fball_sub(r, r, r, n); }
        else fball_sub(r, a, a, n);
        arf_zero(tr);
        check_contains(r, tr, "sub a a", iter);

        /* close operands: b near a (the same exponent, the top limbs
           equal or nearly so), to exercise the order of a difference */
        if (n_randint(state, 4) == 0 && a->size > 0)
        {
            fball_t d;
            fball_init(d);
            fball_set(b, a);
            if (n_randint(state, 2))
            {
                fball_set_ui(d, 1 + n_randint(state, 100));
                fball_mul_2exp_si(d, (a->exp - a->size - n_randint(state, 3)) * FLINT_BITS);
                if (n_randint(state, 2))
                    fball_add(b, b, d, a->size + 4);
                else
                    fball_sub(b, b, d, a->size + 4);
            }
            else if (n_randint(state, 2))
                fball_mul_2exp_si(b, -(slong) FLINT_BITS * (slong) (1 + n_randint(state, 2)));
            if (n_randint(state, 2))
                b->negative = !b->negative;
            fball_random_point(tb, b, state);
            if (aliased == 1) { fball_set(r, a); fball_sub(r, r, b, n); }
            else if (aliased == 2) { fball_set(r, b); fball_sub(r, a, r, n); }
            else fball_sub(r, a, b, n);
            arf_sub(tr, ta, tb, ARF_PREC_EXACT, ARF_RND_NEAR);
            check_contains(r, tr, "sub close", iter);
            if (aliased == 1) { fball_set(r, b); fball_sub(r, r, a, n); }
            else fball_sub(r, b, a, n);
            arf_sub(tr, tb, ta, ARF_PREC_EXACT, ARF_RND_NEAR);
            check_contains(r, tr, "sub close rev", iter);
            fball_clear(d);
        }

        /* mul */
        if (aliased == 1) { fball_set(r, a); fball_mul(r, r, b, n); }
        else if (aliased == 2) { fball_set(r, b); fball_mul(r, a, r, n); }
        else fball_mul(r, a, b, n);
        arf_mul(tr, ta, tb, ARF_PREC_EXACT, ARF_RND_NEAR);
        check_contains(r, tr, "mul", iter);

        /* squaring dispatch */
        fball_mul(r, a, a, n);
        arf_mul(tr, ta, ta, ARF_PREC_EXACT, ARF_RND_NEAR);
        check_contains(r, tr, "sqr", iter);

        /* sqrt / rsqrt / mul_2exp on positive x */
        if (a->size >= 2)
        {
            fball_t xp;
            arb_t rr, tt;
            fball_init(xp);
            arf_t tp;
            arf_init(tp);
            fball_set(xp, a);
            xp->negative = 0;
            if (xp->err > 3)
                xp->err = 3;
            fball_random_point(tp, xp, state);
            if (arf_sgn(tp) > 0)
            {
                arb_init(rr); arb_init(tt);

                fball_sqrt(r, xp, n);
                fball_get_arb(rr, r);
                arb_set_arf(tt, tp);
                arb_sqrt(tt, tt, FLINT_BITS * (n + 6));
                if (!arb_overlaps(rr, tt))
                {
                    flint_printf("FAIL: sqrt (iter %wd)\n", iter);
                    fball_print(xp); fball_print(r);
                    arb_printd(tt, 30); flint_printf("\n");
                    flint_abort();
                }

                fball_rsqrt(r, xp, n);
                fball_get_arb(rr, r);
                arb_set_arf(tt, tp);
                arb_rsqrt(tt, tt, FLINT_BITS * (n + 6));
                if (!arb_overlaps(rr, tt))
                {
                    flint_printf("FAIL: rsqrt (iter %wd)\n", iter);
                    fball_print(xp); fball_print(r);
                    flint_abort();
                }

                arb_clear(rr); arb_clear(tt);
            }
            arf_clear(tp);
            fball_clear(xp);
        }

        /* mul_2exp_si */
        {
            slong sh = (slong) n_randint(state, 300) - 150;
            fball_set(r, a);
            fball_mul_2exp_si(r, sh);
            arf_mul_2exp_si(tr, ta, sh);
            check_contains(r, tr, "mul_2exp", iter);
        }

        /* mul_ui (in place when aliased) */
        {
            ulong c = n_randtest(state);
            if (aliased == 1) { fball_set(r, a); fball_mul_ui(r, r, c, n); }
            else fball_mul_ui(r, a, c, n);
            arf_mul_ui(tr, ta, c, ARF_PREC_EXACT, ARF_RND_NEAR);
            check_contains(r, tr, "mul_ui", iter);
        }

        /* div_ui (in place when aliased) */
        {
            ulong c = n_randtest(state);
            arb_t rr, tt;
            if (c == 0)
                c = 1;
            arb_init(rr); arb_init(tt);
            if (aliased == 1) { fball_set(r, a); fball_div_ui(r, r, c, n); }
            else fball_div_ui(r, a, c, n);
            fball_get_arb(rr, r);
            arb_set_arf(tt, ta);
            arb_div_ui(tt, tt, c, FLINT_BITS * (n + 6));
            if (!arb_overlaps(rr, tt))
            {
                flint_printf("FAIL: div_ui (iter %wd), c = %wu\n", iter, c);
                fball_print(a); fball_print(r);
                flint_abort();
            }
            arb_clear(rr); arb_clear(tt);
        }

        /* div: denominator bounded away from zero, small radius */
        if (b->size >= 2)
        {
            if (b->err != 0)
                b->err = FLINT_MIN(b->err, 3);
            fball_random_point(tb, b, state);
            if (aliased == 1) { fball_set(r, a); fball_div(r, r, b, n); }
            else if (aliased == 2) { fball_set(r, b); fball_div(r, a, r, n); }
            else fball_div(r, a, b, n);
            arf_div(tr, ta, tb, FLINT_BITS * (n + 20), ARF_RND_NEAR);
            /* tr itself is rounded; widen the check with an arb */
            {
                arb_t rr, tt;
                arb_init(rr); arb_init(tt);
                fball_get_arb(rr, r);
                arb_set_arf(tt, tr);
                arb_add_error_2exp_si(tt,
                    arf_is_zero(tr) ? -FLINT_BITS * (n + 19)
                    : (slong) (ARF_EXP(tr) - FLINT_BITS * (n + 19)));
                if (!arb_overlaps(rr, tt))
                {
                    flint_printf("FAIL: div (iter %wd)\n", iter);
                    fball_print(r);
                    arb_printd(tt, 30); flint_printf("\n");
                    flint_abort();
                }
                arb_clear(rr); arb_clear(tt);
            }
        }

        fball_clear(a); fball_clear(b); fball_clear(r);
        arf_clear(ta); arf_clear(tb); arf_clear(tr);
    }
}

static void
test_rsqrt(flint_rand_t state, slong iters)
{
    slong iter;
    for (iter = 0; iter < iters; iter++)
    {
        fball_t r;
        arb_t rr, tt;
        ulong c = 2 + n_randint(state, 1000000);
        slong n = 1 + n_randint(state, 12);

        fball_init(r);
        arb_init(rr); arb_init(tt);

        fball_rsqrt_ui(r, c, n);
        fball_get_arb(rr, r);
        arb_set_ui(tt, c);
        arb_rsqrt(tt, tt, FLINT_BITS * (n + 4));

        if (!arb_overlaps(rr, tt))
        {
            flint_printf("FAIL: rsqrt_ui c=%wu n=%wd\n", c, n);
            fball_print(r);
            flint_abort();
        }
        fball_clear(r);
        arb_clear(rr); arb_clear(tt);
    }
}

static void
test_log2(void)
{
    slong ns[] = { 2, 3, 5, 16, 100, 500 };
    slong i;

    for (i = 0; i < 6; i++)
    {
        slong n = ns[i];
        fball_t v;
        arb_t p1, p2;

        fball_init(v);
        arb_init(p1);
        arb_init(p2);

        fball_const_log2(v, n);
        fball_get_arb(p1, v);
        arb_const_log2(p2, FLINT_BITS * n + 64);

        if (!arb_overlaps(p1, p2))
        {
            flint_printf("FAIL: log2 n=%wd\n", n);
            fball_print(v);
            flint_abort();
        }

        {
            mag_t rad;
            mag_init(rad);
            mag_set_ui(rad, FLINT_MAX(v->err, 1));
            mag_mul_2exp_si(rad, rad,
                FLINT_BITS * (v->exp - v->size));
            if (mag_cmp_2exp_si(rad, -FLINT_BITS * (n - 4)) > 0)
            {
                flint_printf("FAIL: log2 radius too large n=%wd "
                    "err=%wu size=%wd\n", n, v->err, v->size);
                flint_abort();
            }
            mag_clear(rad);
        }

        fball_clear(v);
        arb_clear(p1);
        arb_clear(p2);
    }
}

static void
test_pi(void)
{
    slong ns[] = {2, 3, 4, 5, 8, 16, 33, 100, 331, 1000};
    slong i;

    for (i = 0; i < 10; i++)
    {
        slong n = ns[i];
        fball_t pi;
        arb_t p1, p2;

        fball_init(pi);
        arb_init(p1);
        arb_init(p2);

        fball_const_pi_chudnovsky(pi, n);
        fball_get_arb(p1, pi);
        arb_const_pi(p2, FLINT_BITS * n + 64);

        if (!arb_overlaps(p1, p2))
        {
            flint_printf("FAIL: pi n=%wd\n", n);
            fball_print(pi);
            flint_abort();
        }

        /* the radius should be within a couple of guard limbs */
        {
            mag_t rad;
            mag_init(rad);
            mag_set_ui(rad, FLINT_MAX(pi->err, 1));
            mag_mul_2exp_si(rad, rad,
                FLINT_BITS * (pi->exp - pi->size));
            if (mag_cmp_2exp_si(rad, -FLINT_BITS * (n - 4)) > 0)
            {
                flint_printf("FAIL: pi radius too large n=%wd err=%wu "
                    "size=%wd\n", n, pi->err, pi->size);
                flint_abort();
            }
            mag_clear(rad);
        }

        fball_clear(pi);
        arb_clear(p1);
        arb_clear(p2);
    }
}

/* the small API entry points: setters, swap, exact-zero predicate,
   mpn import, fixed-point export with its rigorous bound, and the
   explicit error-widening helpers */
static void
test_api(flint_rand_t state, slong iters)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        fball_t x, y;
        arb_t ax, ay;
        arf_t t;

        fball_init(x); fball_init(y);
        arb_init(ax); arb_init(ay);
        arf_init(t);

        /* set_ui / set_si / is_zero_exact */
        {
            ulong c = n_randtest(state);
            slong d = (slong) n_randtest(state);

            fball_set_ui(x, c);
            fball_get_arb(ax, x);
            if (!arb_equal_si(ax, (slong) 0) && c == 0)
                { flint_printf("set_ui zero\n"); flint_abort(); }
            arb_set_ui(ay, c);
            if (!arb_equal(ax, ay))
                { flint_printf("set_ui: c = %wu\n", c); flint_abort(); }
            if (!fball_is_zero_exact(x) != !(c == 0))
                { flint_printf("is_zero_exact ui\n"); flint_abort(); }

            fball_set_si(x, d);
            fball_get_arb(ax, x);
            arb_set_si(ay, d);
            if (!arb_equal(ax, ay))
                { flint_printf("set_si: d = %wd\n", d); flint_abort(); }
            fball_add_error_ulps(x, 2.0);
            if (fball_is_zero_exact(x))
                { flint_printf("is_zero_exact with radius\n"); flint_abort(); }
        }

        /* swap */
        {
            fball_randtest(x, state, 6, 1);
            fball_randtest(y, state, 6, 1);
            fball_get_arb(ax, x);
            fball_get_arb(ay, y);
            fball_swap(x, y);
            {
                arb_t bx, by;
                arb_init(bx); arb_init(by);
                fball_get_arb(bx, x);
                fball_get_arb(by, y);
                if (!arb_equal(bx, ay) || !arb_equal(by, ax))
                    { flint_printf("swap\n"); flint_abort(); }
                arb_clear(bx); arb_clear(by);
            }
        }

        /* set_mpn_2exp: exact import */
        {
            ulong p[4];
            slong len = 1 + (slong) n_randint(state, 4);
            slong ebits = (slong) n_randint(state, 400) - 200;
            slong i;

            for (i = 0; i < len; i++)
                p[i] = n_randtest(state);
            fball_set_mpn_2exp(x, p, len, ebits);
            fball_get_arb(ax, x);
            {
                fmpz_t f;
                fmpz_init(f);
                fmpz_set_ui_array(f, p, len);
                arb_set_fmpz(ay, f);
                arb_mul_2exp_si(ay, ay, ebits);
                if (!arb_equal(ax, ay))
                    { flint_printf("set_mpn_2exp: len = %wd, "
                        "ebits = %wd\n", len, ebits); flint_abort(); }
                fmpz_clear(f);
            }
        }

        /* add_error / add_error_ulps: the original point value must
           remain contained after widening */
        {
            fball_randtest(x, state, 6, 1);
            fball_random_point(t, x, state);
            fball_add_error_ulps(x, (double) (1 + n_randint(state, 100)));
            check_contains(x, t, "add_error_ulps", iter);
            fball_add_error(x, (double) (1 + n_randint(state, 100)),
                x->exp - x->size - (slong) n_randint(state, 3));
            check_contains(x, t, "add_error", iter);
        }

        /* get_fixed: build a ball inside [0, 1), export at wn limbs,
           and check the true point against the returned ulp bound */
        {
            ulong p[3], f[6];
            slong len = 1 + (slong) n_randint(state, 3);
            slong wn = 1 + (slong) n_randint(state, 5);
            slong i;
            double bound;

            for (i = 0; i < len; i++)
                p[i] = n_randtest(state);
            p[len - 1] |= UWORD(1) << (FLINT_BITS - 1);
            /* value in [1/2, 1) * 2^-shift */
            fball_set_mpn_2exp(x, p, len,
                -FLINT_BITS * len - (slong) n_randint(state, 40));
            if (n_randint(state, 2))
                fball_add_error_ulps(x,
                    (double) (1 + n_randint(state, 50)));
            fball_random_point(t, x, state);
            if (arf_sgn(t) < 0)
                arf_zero(t);

            bound = fball_get_fixed(f, wn, x);

            {
                arf_t u, w;
                arf_init(u); arf_init(w);
                {
                    fmpz_t g;
                    fmpz_init(g);
                    fmpz_set_ui_array(g, f, wn);
                    arf_set_fmpz(u, g);
                    arf_mul_2exp_si(u, u, -FLINT_BITS * wn);
                    fmpz_clear(g);
                }
                arf_sub(u, t, u, ARF_PREC_EXACT, ARF_RND_NEAR);
                arf_abs(u, u);
                arf_set_d(w, bound);
                arf_mul_2exp_si(w, w, -FLINT_BITS * wn);
                if (arf_cmp(u, w) > 0)
                    { flint_printf("get_fixed bound: wn = %wd, "
                        "bound = %g\n", wn, bound); flint_abort(); }
                arf_clear(u); arf_clear(w);
            }
        }

        /* get_fixed_floor: on success the output must equal the
           EXACT floor of every point of the ball -- checked against
           the floor of a random true point via arf -- and a radius
           spanning a grid line must be rejected */
        {
            ulong p[3], f[6];
            slong len = 1 + (slong) n_randint(state, 3);
            slong wn = 1 + (slong) n_randint(state, 5);
            slong i;

            for (i = 0; i < len; i++)
                p[i] = n_randtest(state);
            p[len - 1] |= UWORD(1) << (FLINT_BITS - 1);
            fball_set_mpn_2exp(x, p, len,
                -FLINT_BITS * len - (slong) n_randint(state, 3));
            if (n_randint(state, 2))
                fball_add_error_ulps(x,
                    (double) (1 + n_randint(state, 1000)));

            if (fball_get_fixed_floor(f, wn, x))
            {
                fmpz_t g, h;
                fmpz_init(g); fmpz_init(h);
                fball_random_point(t, x, state);
                if (arf_sgn(t) < 0)
                    arf_zero(t);
                arf_mul_2exp_si(t, t, FLINT_BITS * wn);
                arf_get_fmpz(g, t, ARF_RND_FLOOR);
                fmpz_set_ui_array(h, f, wn);
                if (!fmpz_equal(g, h))
                    { flint_printf("FAIL: get_fixed_floor value "
                        "(iter %wd)\n", iter); flint_abort(); }
                fmpz_clear(g); fmpz_clear(h);
            }

            /* a ball straddling a grid line must be refused: center
               the value ON a grid multiple with a nonzero radius */
            fball_set_mpn_2exp(x, p, 1, -FLINT_BITS);
            fball_add_error_ulps(x, 1.0);
            if (p[0] != 0 && fball_get_fixed_floor(f, 1, x))
                { flint_printf("FAIL: get_fixed_floor accepted a "
                    "grid-straddling ball (iter %wd)\n", iter);
                  flint_abort(); }
        }

        fball_clear(x); fball_clear(y);
        arb_clear(ax); arb_clear(ay);
        arf_clear(t);
    }
}

/* fball_submul_bounded against fball_mul + fball_sub: a = b c + d with
   |d| < B^E, at various precisions, exact and inexact operands */
static void
test_submul(flint_rand_t state, slong iters)
{
    slong iter;
    for (iter = 0; iter < iters; iter++)
    {
        fball_t a, b, c, d, p, r1, r2;
        arb_t x, y, xa;
        slong n = 1 + n_randint(state, (iter % 10 == 0) ? 300 : 20);
        slong sb = 1 + n_randint(state, 2 * n + 1), sc = 1 + n_randint(state, 2 * n + 1);
        slong E, sd;
        nn_ptr tmp;

        fball_init(a); fball_init(b); fball_init(c); fball_init(d);
        fball_init(p); fball_init(r1); fball_init(r2);
        arb_init(x); arb_init(y); arb_init(xa);

        tmp = flint_malloc((2 * n + 2) * sizeof(ulong));
        flint_mpn_rrandom(tmp, state, sb);
        if (tmp[sb - 1] == 0) tmp[sb - 1] = 1;
        fball_set_mpn_2exp(b, tmp, sb, (slong) n_randint(state, 300) - 150);
        flint_mpn_rrandom(tmp, state, sc);
        if (tmp[sc - 1] == 0) tmp[sc - 1] = 1;
        fball_set_mpn_2exp(c, tmp, sc, (slong) n_randint(state, 300) - 150);
        if (n_randint(state, 2)) fball_neg(b);
        if (n_randint(state, 2)) fball_neg(c);
        if (n_randint(state, 3) == 0) fball_add_error_ulps(b, 1.0 + n_randint(state, 5));
        if (n_randint(state, 3) == 0) fball_add_error_ulps(c, 1.0 + n_randint(state, 5));

        /* the exact-ish product, and a small d */
        fball_mul(p, b, c, sb + sc + 2);
        sd = 1 + n_randint(state, n + 1);
        flint_mpn_rrandom(tmp, state, sd);
        if (tmp[sd - 1] == 0) tmp[sd - 1] = 1;
        /* d = (tmp, sd) B^(p.exp - n - 1 - sd + t) < B^(p.exp - n - 1 + t) */
        {
            slong t = n_randint(state, 3);
            fball_set_mpn_2exp(d, tmp, sd, FLINT_BITS * (p->exp - n - 1 - sd + t));
            E = p->exp - n - 1 + t;
        }
        if (n_randint(state, 2)) fball_neg(d);
        if (n_randint(state, 5) == 0) fball_zero(d);
        /* a = p + d, exact when p is */
        fball_add(a, p, d, sb + sc + n + 4);
        if (n_randint(state, 3) == 0) fball_add_error_ulps(a, 1.0 + n_randint(state, 5));
        /* the bound |a - b c| < B^E: d's magnitude, plus the radii,
           which are far below it */
        E = E + 1;

        fball_submul_bounded(r1, a, b, c, E, n);
        fball_mul(r2, b, c, sb + sc + 2);
        fball_sub(r2, a, r2, sb + sc + n + 4);

        fball_get_arb(x, r1);
        fball_get_arb(y, r2);
        if (!arb_overlaps(x, y))
        {
            flint_printf("FAIL: submul_bounded (iter %wd), n = %wd, E = %wd\n",
                iter, n, E);
            fball_print(a); fball_print(b); fball_print(c);
            fball_print(r1); fball_print(r2);
            flint_abort();
        }
        /* accuracy: the limbs of the residual above B^(E - n), its
           weight window, less a partial top limb and the 3 units */
        fball_get_arb(xa, a);
        if (a->err == 0 && b->err == 0 && c->err == 0
            && r1->size > 0
            && arb_rel_accuracy_bits(x) < FLINT_BITS * (r1->exp - (E - n) - 1) - 8)
        {
            flint_printf("FAIL: submul_bounded accuracy (iter %wd), n = %wd: "
                "%wd bits\n", iter, n, arb_rel_accuracy_bits(x));
            fball_print(r1); fball_print(r2);
            flint_abort();
        }

        flint_free(tmp);
        fball_clear(a); fball_clear(b); fball_clear(c); fball_clear(d);
        fball_clear(p); fball_clear(r1); fball_clear(r2);
        arb_clear(x); arb_clear(y); arb_clear(xa);
    }
}

/* the complex product against the exact arf products, with the parts
   of an operand at unrelated magnitudes and every aliasing */
static void
test_mul_complex(flint_rand_t state, slong iters)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        fball_t ar, ai, br, bi, rr, ri, t1, t2;
        arf_t xr, xi, yr, yi, zr, zi, u;
        slong n = 1 + n_randint(state, (iter % 10 == 0) ? 400 : 12);
        slong maxsize = (iter % 10 == 0) ? 400 : 10;
        int alias = (int) n_randint(state, 5);

        fball_init(ar); fball_init(ai); fball_init(br); fball_init(bi);
        fball_init(rr); fball_init(ri); fball_init(t1); fball_init(t2);
        arf_init(xr); arf_init(xi); arf_init(yr); arf_init(yi);
        arf_init(zr); arf_init(zi); arf_init(u);

        fball_randtest(ar, state, maxsize, 1);
        fball_randtest(ai, state, maxsize, 1);
        fball_randtest(br, state, maxsize, 1);
        fball_randtest(bi, state, maxsize, 1);
        /* often a tiny or an exact-zero imaginary part, as in the
           accumulation of exp(i x) for small x */
        if (n_randint(state, 3) == 0)
            fball_mul_2exp_si(ai, -(slong) FLINT_BITS * (slong) n_randint(state, 2 * maxsize + 2));
        if (n_randint(state, 5) == 0)
            fball_zero(ai);
        if (n_randint(state, 3) == 0)
            fball_mul_2exp_si(bi, -(slong) FLINT_BITS * (slong) n_randint(state, 2 * maxsize + 2));
        if (n_randint(state, 5) == 0)
            fball_zero(bi);
        fball_random_point(xr, ar, state);
        fball_random_point(xi, ai, state);
        fball_random_point(yr, br, state);
        fball_random_point(yi, bi, state);

        arf_mul(zr, xr, yr, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul(u, xi, yi, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_sub(zr, zr, u, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul(zi, xr, yi, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul(u, xi, yr, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_add(zi, zi, u, ARF_PREC_EXACT, ARF_RND_NEAR);

        if (alias == 0)
            fball_mul_complex(rr, ri, ar, ai, br, bi, n);
        else if (alias == 1)
        {
            fball_set(rr, ar); fball_set(ri, ai);
            fball_mul_complex(rr, ri, rr, ri, br, bi, n);
        }
        else if (alias == 2)
        {
            fball_set(rr, br); fball_set(ri, bi);
            fball_mul_complex(rr, ri, ar, ai, rr, ri, n);
        }
        else if (alias == 3)
        {
            /* a square */
            fball_set(br, ar); fball_set(bi, ai);
            arf_set(yr, xr); arf_set(yi, xi);
            arf_mul(zr, xr, yr, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_mul(u, xi, yi, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_sub(zr, zr, u, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_mul(zi, xr, yi, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_mul(u, xi, yr, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_add(zi, zi, u, ARF_PREC_EXACT, ARF_RND_NEAR);
            fball_mul_complex(rr, ri, ar, ai, ar, ai, n);
        }
        else
        {
            fball_set(rr, ai); fball_set(ri, ar);
            fball_mul_complex(rr, ri, ri, rr, br, bi, n);
        }

        check_contains(rr, zr, "mul_complex re", iter);
        check_contains(ri, zi, "mul_complex im", iter);

        /* accuracy: about n limbs relative to the larger part of the
           result (one frame for both parts), or that of the operands:
           against the four products, the radius of each part is
           within a few limbs of the reference's or of the frame */
        {
            fball_t t3, t4;
            double lr, li, ref_r, ref_i, top, fr;
            fball_init(t3); fball_init(t4);
            fball_mul(t1, ar, br, n);
            fball_mul(t2, ai, bi, n);
            fball_sub(t1, t1, t2, n);
            fball_mul(t3, ar, bi, n);
            fball_mul(t4, ai, br, n);
            fball_add(t3, t3, t4, n);
            /* log2 of a radius, and of the frame's ulp */
#define RADBITS(x) ((x)->err == 0 ? -1e300 : (double) FLINT_BITS * ((x)->exp - (x)->size) + FLINT_BIT_COUNT((x)->err))
#define TOPBITS(x) ((x)->size == 0 ? -1e300 : (double) FLINT_BITS * (x)->exp)
            lr = RADBITS(rr); li = RADBITS(ri);
            ref_r = RADBITS(t1); ref_i = RADBITS(t3);
            top = FLINT_MAX(TOPBITS(t1), TOPBITS(t3));
            fr = top - FLINT_BITS * (n - 2);
            if (lr > FLINT_MAX(FLINT_MAX(ref_r, ref_i), fr) + 3 * FLINT_BITS
                || li > FLINT_MAX(FLINT_MAX(ref_r, ref_i), fr) + 3 * FLINT_BITS)
            {
                flint_printf("FAIL: mul_complex accuracy (iter %wd): rad %g %g ref %g %g frame %g\n", iter, lr, li, ref_r, ref_i, fr);
                fball_print(ar); fball_print(ai); fball_print(br); fball_print(bi);
                fball_print(rr); fball_print(ri); fball_print(t1); fball_print(t3);
                flint_abort();
            }
            fball_clear(t3); fball_clear(t4);
        }

        fball_clear(ar); fball_clear(ai); fball_clear(br); fball_clear(bi);
        fball_clear(rr); fball_clear(ri); fball_clear(t1); fball_clear(t2);
        arf_clear(xr); arf_clear(xi); arf_clear(yr); arf_clear(yi);
        arf_clear(zr); arf_clear(zi); arf_clear(u);
    }
}

/* the AGM against arb_agm at the midpoints, for exact and inexact
   inputs, equal inputs, nearly equal inputs, far-apart inputs and
   every order of the finish */
static void
test_agm(flint_rand_t state, slong iters)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        slong n = 1 + n_randint(state, (iter % 10 == 0) ? 300 : 30);
        int m = n_randint(state, 3) ? 0 : 2 + (int) n_randint(state, 15);
        fball_t x, y, r;
        arb_t xa, ya, ra, rr;
        int w;

        fball_init(x); fball_init(y); fball_init(r);
        arb_init(xa); arb_init(ya); arb_init(ra); arb_init(rr);

        for (w = 0; w < 2; w++)
        {
            fball_struct * v = w ? y : x;
            slong sz = 1 + n_randint(state, n + 2);
            nn_ptr t = flint_malloc(sz * sizeof(ulong));
            flint_mpn_urandomb(t, state, FLINT_BITS * sz);
            t[sz - 1] |= 1;
            if (n_randint(state, 5) == 0)
                fball_set_ui(v, 1 + n_randint(state, 10));
            else
                fball_set_mpn_2exp(v, t, sz, -FLINT_BITS * sz
                    + (slong) n_randint(state, 400) - 200);
            if (n_randint(state, 3) == 0)
                fball_add_error_ulps(v, 1 + n_randint(state, 100));
            flint_free(t);
        }
        if (n_randint(state, 6) == 0)
            fball_set(y, x);
        if (n_randint(state, 6) == 0)
        {
            fball_set(y, x);
            fball_add_error_2exp(y, fball_mag_2exp(x)
                - FLINT_BITS * (slong) n_randint(state, n + 3));
        }

        _fball_agm_order(r, x, y, n, m);
        fball_get_arb(xa, x);
        fball_get_arb(ya, y);
        fball_get_arb(rr, r);
        arb_get_mid_arb(xa, xa);
        arb_get_mid_arb(ya, ya);
        arb_agm(ra, xa, ya, FLINT_BITS * n + 128);
        if (!arb_overlaps(ra, rr))
        {
            flint_printf("FAIL: agm (iter %wd, n = %wd, m = %d)\n", iter, n, m);
            fball_print(x); fball_print(y); fball_print(r);
            flint_abort();
        }
        if (x->err == 0 && y->err == 0
            && arb_rel_accuracy_bits(rr) < FLINT_BITS * n - 16)
        {
            flint_printf("FAIL: agm accuracy (iter %wd, n = %wd, m = %d): %wd bits\n",
                iter, n, m, arb_rel_accuracy_bits(rr));
            fball_print(x); fball_print(y); fball_print(r);
            flint_abort();
        }

        fball_clear(x); fball_clear(y); fball_clear(r);
        arb_clear(xa); arb_clear(ya); arb_clear(ra); arb_clear(rr);
    }
}

TEST_FUNCTION_START(fixed_fball, state)
{
    test_ops(state, 3000 + 3000 * flint_test_multiplier(), 10, 12);
    /* larger operands: the middle-product path of fball_mul and the
       divide-and-conquer divisions */
    test_ops(state, 200 + 200 * flint_test_multiplier(), 120, 100);
    /* a few operands long enough for the Newton divisions and square
       roots (FBALL_DIV_NEWTON_CUTOFF, FBALL_SQRT_NEWTON_CUTOFF) */
    test_ops(state, 12 * flint_test_multiplier(), 2400, 2400);
    test_api(state, 1000 + 1000 * flint_test_multiplier());
    test_rsqrt(state, 200 + 200 * flint_test_multiplier());
    test_submul(state, 300 + 300 * flint_test_multiplier());
    test_mul_complex(state, 2000 + 2000 * flint_test_multiplier());
    test_agm(state, 300 + 300 * flint_test_multiplier());
    test_pi();
    test_log2();

    TEST_FUNCTION_END(state);
}
