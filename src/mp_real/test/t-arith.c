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
#include "mp_real.h"
#include "mp_real/impl.h"
#include "ball_helpers.h"

/* the arithmetic on random balls: every output ball must contain the
   exact result of a random true point drawn inside each operand ball
   (with aliasing variants), and the normalization invariant must hold
   afterward */

static void
test_ops(flint_rand_t state, slong iters, slong maxsize, slong maxn)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        mp_real_t a, b, r;
        arf_t ta, tb, tr;
        slong n = 1 + n_randint(state, maxn);
        int aliased = (int) n_randint(state, 3);

        mp_real_init(a); mp_real_init(b); mp_real_init(r);
        arf_init(ta); arf_init(tb); arf_init(tr);

        mp_real_randtest(a, state, maxsize, 1);
        mp_real_randtest(b, state, maxsize, 1);
        mp_real_random_point(ta, a, state);
        mp_real_random_point(tb, b, state);

        /* add */
        if (aliased == 1) { mp_real_set(r, a); mp_real_add(r, r, b, n); }
        else if (aliased == 2) { mp_real_set(r, b); mp_real_add(r, a, r, n); }
        else mp_real_add(r, a, b, n);
        arf_add(tr, ta, tb, ARF_PREC_EXACT, ARF_RND_NEAR);
        check_contains(r, tr, "add", iter);

        /* sub (aliased variants exercise the in-place window
           build and negation) */
        if (aliased == 1) { mp_real_set(r, a); mp_real_sub(r, r, b, n); }
        else if (aliased == 2) { mp_real_set(r, b); mp_real_sub(r, a, r, n); }
        else mp_real_sub(r, a, b, n);
        arf_sub(tr, ta, tb, ARF_PREC_EXACT, ARF_RND_NEAR);
        check_contains(r, tr, "sub", iter);

        /* addmul_ui / submul_ui, all aliasings; b is often placed
           inside a's window (the in-place path) */
        {
            ulong c = n_randint(state, 4) ? n_randtest(state) : n_randint(state, 10);
            mp_real_t bb;
            arf_t tbb, u;
            int sub = n_randint(state, 2);
            mp_real_init(bb); arf_init(tbb); arf_init(u);
            mp_real_set(bb, b);
            if (n_randint(state, 2) && a->size > 2 && b->size > 0)
            {
                /* shift b to sit within a's limbs, or to reach below
                   a's bottom while ending under its top */
                mp_real_mul_2exp_si(bb, bb, FLINT_BITS * ((a->exp - 2 - (slong) n_randint(state, 2)) - bb->exp));
                if (n_randint(state, 2))
                    mp_real_mul_2exp_si(bb, bb, -FLINT_BITS * (slong) n_randint(state, 3));
            }
            mp_real_random_point(tbb, bb, state);
            arf_mul_ui(u, tbb, c, ARF_PREC_EXACT, ARF_RND_NEAR);
            if (sub)
                arf_sub(tr, ta, u, ARF_PREC_EXACT, ARF_RND_NEAR);
            else
                arf_add(tr, ta, u, ARF_PREC_EXACT, ARF_RND_NEAR);
            if (aliased == 1) { mp_real_set(r, a); if (sub) mp_real_submul_ui(r, r, bb, c, n); else mp_real_addmul_ui(r, r, bb, c, n); }
            else if (aliased == 2) { mp_real_set(r, bb); if (sub) mp_real_submul_ui(r, a, r, c, n); else mp_real_addmul_ui(r, a, r, c, n); }
            else { if (sub) mp_real_submul_ui(r, a, bb, c, n); else mp_real_addmul_ui(r, a, bb, c, n); }
            check_contains(r, tr, sub ? "submul_ui" : "addmul_ui", iter);
            /* exact inputs fitting the precision: exact result */
            if (a->err == 0 && bb->err == 0 && r->err != 0 && a->size + 2 < n
                && bb->exp - bb->size >= a->exp - a->size && bb->exp + 1 <= a->exp)
            {
                flint_printf("FAIL: addmul_ui inexact (iter %wd)\n", iter);
                mp_real_print(a); mp_real_print(bb); mp_real_print(r);
                flint_abort();
            }
            mp_real_clear(bb); arf_clear(tbb); arf_clear(u);
        }

        /* the same ball twice (doubling and cancellation), also in
           place */
        if (aliased == 1) { mp_real_set(r, a); mp_real_add(r, r, r, n); }
        else mp_real_add(r, a, a, n);
        arf_add(tr, ta, ta, ARF_PREC_EXACT, ARF_RND_NEAR);
        check_contains(r, tr, "add a a", iter);
        if (aliased == 1) { mp_real_set(r, a); mp_real_sub(r, r, r, n); }
        else mp_real_sub(r, a, a, n);
        arf_zero(tr);
        check_contains(r, tr, "sub a a", iter);

        /* close operands: b near a (the same exponent, the top limbs
           equal or nearly so), to exercise the order of a difference */
        if (n_randint(state, 4) == 0 && a->size > 0)
        {
            mp_real_t d;
            mp_real_init(d);
            mp_real_set(b, a);
            if (n_randint(state, 2))
            {
                mp_real_set_ui(d, 1 + n_randint(state, 100));
                mp_real_mul_2exp_si(d, d, (a->exp - a->size - n_randint(state, 3)) * FLINT_BITS);
                if (n_randint(state, 2))
                    mp_real_add(b, b, d, a->size + 4);
                else
                    mp_real_sub(b, b, d, a->size + 4);
            }
            else if (n_randint(state, 2))
                mp_real_mul_2exp_si(b, b, -(slong) FLINT_BITS * (slong) (1 + n_randint(state, 2)));
            if (n_randint(state, 2))
                b->negative = !b->negative;
            mp_real_random_point(tb, b, state);
            if (aliased == 1) { mp_real_set(r, a); mp_real_sub(r, r, b, n); }
            else if (aliased == 2) { mp_real_set(r, b); mp_real_sub(r, a, r, n); }
            else mp_real_sub(r, a, b, n);
            arf_sub(tr, ta, tb, ARF_PREC_EXACT, ARF_RND_NEAR);
            check_contains(r, tr, "sub close", iter);
            if (aliased == 1) { mp_real_set(r, b); mp_real_sub(r, r, a, n); }
            else mp_real_sub(r, b, a, n);
            arf_sub(tr, tb, ta, ARF_PREC_EXACT, ARF_RND_NEAR);
            check_contains(r, tr, "sub close rev", iter);
            mp_real_clear(d);
        }

        /* mul */
        if (aliased == 1) { mp_real_set(r, a); mp_real_mul(r, r, b, n); }
        else if (aliased == 2) { mp_real_set(r, b); mp_real_mul(r, a, r, n); }
        else mp_real_mul(r, a, b, n);
        arf_mul(tr, ta, tb, ARF_PREC_EXACT, ARF_RND_NEAR);
        check_contains(r, tr, "mul", iter);

        /* squaring dispatch */
        mp_real_mul(r, a, a, n);
        arf_mul(tr, ta, ta, ARF_PREC_EXACT, ARF_RND_NEAR);
        check_contains(r, tr, "sqr", iter);

        /* sqrt / rsqrt / mul_2exp on positive x */
        if (a->size >= 2)
        {
            mp_real_t xp;
            arb_t rr, tt;
            mp_real_init(xp);
            arf_t tp;
            arf_init(tp);
            mp_real_set(xp, a);
            xp->negative = 0;
            if (xp->err > 3)
                xp->err = 3;
            mp_real_random_point(tp, xp, state);
            if (arf_sgn(tp) > 0)
            {
                arb_init(rr); arb_init(tt);

                mp_real_sqrt(r, xp, n);
                mp_real_get_arb(rr, r);
                arb_set_arf(tt, tp);
                arb_sqrt(tt, tt, FLINT_BITS * (n + 6));
                if (!arb_overlaps(rr, tt))
                {
                    flint_printf("FAIL: sqrt (iter %wd)\n", iter);
                    mp_real_print(xp); mp_real_print(r);
                    arb_printd(tt, 30); flint_printf("\n");
                    flint_abort();
                }

                mp_real_rsqrt(r, xp, n);
                mp_real_get_arb(rr, r);
                arb_set_arf(tt, tp);
                arb_rsqrt(tt, tt, FLINT_BITS * (n + 6));
                if (!arb_overlaps(rr, tt))
                {
                    flint_printf("FAIL: rsqrt (iter %wd)\n", iter);
                    mp_real_print(xp); mp_real_print(r);
                    flint_abort();
                }

                arb_clear(rr); arb_clear(tt);
            }
            arf_clear(tp);
            mp_real_clear(xp);
        }

        /* mul_2exp_si */
        {
            slong sh = (slong) n_randint(state, 300) - 150;
            mp_real_set(r, a);
            mp_real_mul_2exp_si(r, r, sh);
            arf_mul_2exp_si(tr, ta, sh);
            check_contains(r, tr, "mul_2exp", iter);
        }

        /* mul_ui (in place when aliased) */
        {
            ulong c = n_randtest(state);
            if (aliased == 1) { mp_real_set(r, a); mp_real_mul_ui(r, r, c, n); }
            else mp_real_mul_ui(r, a, c, n);
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
            if (aliased == 1) { mp_real_set(r, a); mp_real_div_ui(r, r, c, n); }
            else mp_real_div_ui(r, a, c, n);
            mp_real_get_arb(rr, r);
            arb_set_arf(tt, ta);
            arb_div_ui(tt, tt, c, FLINT_BITS * (n + 6));
            if (!arb_overlaps(rr, tt))
            {
                flint_printf("FAIL: div_ui (iter %wd), c = %wu\n", iter, c);
                mp_real_print(a); mp_real_print(r);
                flint_abort();
            }
            arb_clear(rr); arb_clear(tt);
        }

        /* div: denominator bounded away from zero, small radius */
        if (b->size >= 2)
        {
            if (b->err != 0)
                b->err = FLINT_MIN(b->err, 3);
            mp_real_random_point(tb, b, state);
            if (aliased == 1) { mp_real_set(r, a); mp_real_div(r, r, b, n); }
            else if (aliased == 2) { mp_real_set(r, b); mp_real_div(r, a, r, n); }
            else mp_real_div(r, a, b, n);
            arf_div(tr, ta, tb, FLINT_BITS * (n + 20), ARF_RND_NEAR);
            /* tr itself is rounded; widen the check with an arb */
            {
                arb_t rr, tt;
                arb_init(rr); arb_init(tt);
                mp_real_get_arb(rr, r);
                arb_set_arf(tt, tr);
                arb_add_error_2exp_si(tt,
                    arf_is_zero(tr) ? -FLINT_BITS * (n + 19)
                    : (slong) (ARF_EXP(tr) - FLINT_BITS * (n + 19)));
                if (!arb_overlaps(rr, tt))
                {
                    flint_printf("FAIL: div (iter %wd)\n", iter);
                    mp_real_print(r);
                    arb_printd(tt, 30); flint_printf("\n");
                    flint_abort();
                }
                arb_clear(rr); arb_clear(tt);
            }
        }

        mp_real_clear(a); mp_real_clear(b); mp_real_clear(r);
        arf_clear(ta); arf_clear(tb); arf_clear(tr);
    }
}

static void
test_rsqrt(flint_rand_t state, slong iters)
{
    slong iter;
    for (iter = 0; iter < iters; iter++)
    {
        mp_real_t r;
        arb_t rr, tt;
        ulong c = 2 + n_randint(state, 1000000);
        slong n = 1 + n_randint(state, 12);

        mp_real_init(r);
        arb_init(rr); arb_init(tt);

        mp_real_rsqrt_ui(r, c, n);
        mp_real_get_arb(rr, r);
        arb_set_ui(tt, c);
        arb_rsqrt(tt, tt, FLINT_BITS * (n + 4));

        if (!arb_overlaps(rr, tt))
        {
            flint_printf("FAIL: rsqrt_ui c=%wu n=%wd\n", c, n);
            mp_real_print(r);
            flint_abort();
        }
        mp_real_clear(r);
        arb_clear(rr); arb_clear(tt);
    }
}

TEST_FUNCTION_START(mp_real_arith, state)
{
    test_ops(state, 3000 + 3000 * flint_test_multiplier(), 10, 12);
    /* larger operands: the middle-product path of mp_real_mul and the
       divide-and-conquer divisions */
    test_ops(state, 200 + 200 * flint_test_multiplier(), 120, 100);
    /* a few operands long enough for the Newton divisions and square
       roots (MP_REAL_BALL_DIV_NEWTON_CUTOFF, MP_REAL_BALL_SQRT_NEWTON_CUTOFF) */
    test_ops(state, 12 * flint_test_multiplier(), 2400, 2400);
    test_rsqrt(state, 200 + 200 * flint_test_multiplier());

    TEST_FUNCTION_END(state);
}
