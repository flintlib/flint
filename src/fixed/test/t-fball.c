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

/* random normalized fball; if allow_err, a random radius */
static void
fball_randtest(fball_t x, flint_rand_t state, slong maxsize, int allow_err)
{
    slong size = n_randint(state, maxsize + 1);

    if (size == 0 && n_randint(state, 2))
    {
        fball_zero(x);
        x->exp = (slong) n_randint(state, 13) - 6;
        if (allow_err && n_randint(state, 2))
            x->err = ldexp((double) (1 + n_randint(state, 1000)),
                (int) n_randint(state, 60) - 30);
        x->erra = 0;
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
    x->err = 0.0;
    x->erra = 0;
    if (allow_err && n_randint(state, 2))
        x->err = ldexp((double) (1 + n_randint(state, 1000)),
            (int) n_randint(state, 60) - 30);
    if (x->err == 0.0 && n_randint(state, 2))
    {
        /* exercise the low-zero-limb stripping */
        x->d[0] &= ~(ulong) n_randint(state, 2);
    }
    if (x->err > FBALL_ERR_MAX)
        x->err = FBALL_ERR_MAX;
    /* keep the normalization invariants */
    if (x->err == 0.0)
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
    if (x->err != 0.0)
    {
        /* mid + (r / 2^30) * err * ulp, r in [-2^30, 2^30] */
        slong r = (slong) n_randint(state, UWORD(1) << 31)
                    - (slong) (UWORD(1) << 30);
        arf_set_d(u, x->err);
        arf_mul_si(u, u, r, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul_2exp_si(u, u,
            FLINT_BITS * (x->exp - x->size + x->erra) - 30);
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
    /* err must remain normalized (unless the mantissa collapsed) */
    if (res->size > 0 && res->err > FBALL_ERR_MAX * 68.0)
    {
        flint_printf("FAIL: %s unnormalized err %g (iter %wd)\n",
            op, res->err, iter);
        flint_abort();
    }
    arb_clear(r);
}

static void
test_ops(flint_rand_t state, slong iters)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        fball_t a, b, r;
        arf_t ta, tb, tr;
        slong n = 1 + n_randint(state, 12);
        int aliased = (int) n_randint(state, 3);

        fball_init(a); fball_init(b); fball_init(r);
        arf_init(ta); arf_init(tb); arf_init(tr);

        fball_randtest(a, state, 10, 1);
        fball_randtest(b, state, 10, 1);
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
        if (a->size >= 2 && (a->err == 0.0 || a->size >= 2))
        {
            fball_t xp;
            arb_t rr, tt;
            fball_init(xp);
            arf_t tp;
            arf_init(tp);
            fball_set(xp, a);
            xp->negative = 0;
            if (xp->err > 3.0)
                xp->err = 3.0;
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

        /* div: denominator bounded away from zero, small radius */
        if (b->size >= 2)
        {
            if (b->err != 0.0)
                b->err = FLINT_MIN(b->err, 3.0);
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
            mag_set_d(rad, FLINT_MAX(v->err, 1.0));
            mag_mul_2exp_si(rad, rad,
                FLINT_BITS * (v->exp - v->size + v->erra));
            if (mag_cmp_2exp_si(rad, -FLINT_BITS * (n - 4)) > 0)
            {
                flint_printf("FAIL: log2 radius too large n=%wd "
                    "err=%g size=%wd\n", n, v->err, v->size);
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
            mag_set_d(rad, FLINT_MAX(pi->err, 1.0));
            mag_mul_2exp_si(rad, rad,
                FLINT_BITS * (pi->exp - pi->size + pi->erra));
            if (mag_cmp_2exp_si(rad, -FLINT_BITS * (n - 4)) > 0)
            {
                flint_printf("FAIL: pi radius too large n=%wd err=%g "
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

TEST_FUNCTION_START(fixed_fball, state)
{
    test_ops(state, 3000 + 3000 * flint_test_multiplier());
    test_api(state, 1000 + 1000 * flint_test_multiplier());
    test_rsqrt(state, 200 + 200 * flint_test_multiplier());
    test_pi();
    test_log2();

    TEST_FUNCTION_END(state);
}
