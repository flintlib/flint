/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Shared checks for the elementary functions of arf: the result against
   a ball enclosing the value (correctly rounded in the strict modes, within
   the tolerance of the relaxed modes), with random arguments in the
   tiny, normal and huge ranges. */

#ifndef ARF_TEST_ELEM_CHECK_H
#define ARF_TEST_ELEM_CHECK_H

#include "test_helpers.h"
#include "arf.h"
#include "arb.h"
#include "gr.h"

typedef int (* arf_elem_func_t)(arf_ptr, arf_srcptr, slong, arf_rnd_t);
typedef void (* arb_elem_func_t)(arb_ptr, const arb_t, slong);

/* a random argument: the mantissa with random bits, the exponent in the
   normal range, or tiny (below -prec), or huge, or a bignum */
static void
arf_elem_randtest_arg(arf_t x, flint_rand_t state, slong prec)
{
    slong bits = 1 + n_randint(state, 300);
    slong e;

    if (n_randint(state, 8) == 0)
        bits = 1 + n_randint(state, 3000);

    if (n_randint(state, 30) == 0)
    {
        arf_randtest_special(x, state, bits, 10);
        return;
    }

    arf_urandom(x, state, bits, ARF_RND_DOWN);
    if (arf_is_zero(x))
        arf_one(x);
    if (n_randint(state, 2))
        arf_neg(x, x);

    switch (n_randint(state, 10))
    {
        case 8:
            /* a power of two, or a few bits, at a wide exponent */
            arf_set_ui(x, 1 + n_randint(state, 4) * n_randint(state, 2));
            arf_mul_2exp_si(x, x, (slong) n_randint(state, 4 * prec + 200) - 2 * prec - 100);
            if (n_randint(state, 2))
                arf_neg(x, x);
            break;
        case 9:
            /* 1 +- m 2^-k, close to one */
            arf_set_ui(x, 1 + n_randint(state, 8));
            arf_mul_2exp_si(x, x, -(slong) (1 + n_randint(state, 2 * prec + 40)));
            if (n_randint(state, 2))
                arf_neg(x, x);
            arf_add_ui(x, x, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
            if (n_randint(state, 2))
                arf_neg(x, x);
            break;
        case 0:
            /* tiny */
            e = -prec - (slong) n_randint(state, 100);
            arf_mul_2exp_si(x, x, e);
            break;
        case 1:
            /* around the tiny threshold */
            e = -prec + 20 - (slong) n_randint(state, 40);
            arf_mul_2exp_si(x, x, e);
            break;
        case 2:
            /* huge */
            e = 50 + n_randint(state, 100 + 2 * prec);
            arf_mul_2exp_si(x, x, e);
            break;
        case 3:
            {
                /* a bignum exponent */
                fmpz_t t;
                fmpz_init(t);
                fmpz_randtest(t, state, 200);
                fmpz_mul_2exp(t, t, 60);
                arf_mul_2exp_fmpz(x, x, t);
                fmpz_clear(t);
            }
            break;
        default:
            e = (slong) n_randint(state, 200) - 100;
            arf_mul_2exp_si(x, x, e);
            break;
    }
}

/* Checks the result r = f(x) (with the inexact flag) at prec bits in the
   mode rnd against the ball b enclosing the value, computed at prec + 64
   bits: in the strict modes the rounding of the ball's endpoints must
   agree with r (if the two roundings differ the check is inconclusive),
   and the inexact flag must be as given; in the relaxed modes
   |r - mid| <= rad + tol ulps of r.  A nan is allowed for a non-finite
   x, a nan reference or |x| >= 2^maglim. */
static void
arf_elem_verify2(const char * name, const arf_t x, const arf_t y, const arf_t r, int inexact,
    const arb_t b, slong prec, arf_rnd_t rnd, slong maglim)
{
    arf_t s, lo, hi, u;
    int ok, exact_expected, exact_got;

    arf_init(s);
    arf_init(lo);
    arf_init(hi);
    arf_init(u);

    if (arf_is_nan(r))
    {
        /* allowed only beyond the magnitude limit, for a nan argument or
           outside the domain (where arb is indeterminate) */
        if (!(!arf_is_finite(x) || arf_is_nan(arb_midref(b))
                || arf_cmpabs_2exp_si(x, maglim) >= 0))
        {
            flint_printf("FAIL (nan): %s, prec = %wd, rnd = %d\n", name, prec, rnd);
            flint_printf("x = "); arf_debug(x); flint_printf("\n");
        if (y != NULL) { flint_printf("y = "); arf_debug(y); flint_printf("\n"); }
            flint_printf("b = "); arb_printd(b, 30); flint_printf("\n");
            flint_abort();
        }
        goto cleanup;
    }

    if (!arb_is_finite(b))
        goto cleanup;

    exact_got = (inexact == 0);
    exact_expected = arb_is_exact(b);

    if (!arf_is_finite(r))
    {
        /* an infinite result for an infinite argument only */
        if (!arb_is_exact(b) || !arf_equal(arb_midref(b), r))
        {
            flint_printf("FAIL (infinite): %s, prec = %wd, rnd = %d\n", name, prec, rnd);
            flint_printf("x = "); arf_debug(x); flint_printf("\n");
        if (y != NULL) { flint_printf("y = "); arf_debug(y); flint_printf("\n"); }
            flint_printf("r = "); arf_debug(r); flint_printf("\n");
            flint_abort();
        }
        goto cleanup;
    }

    arb_get_interval_arf(lo, hi, b, prec + 200);

    if (rnd < ARF_RND_FAST)
    {
        arf_set_round(lo, lo, prec, rnd);
        arf_set_round(hi, hi, prec, rnd);

        /* (an undecided reference, whose endpoints round differently,
           only bounds the result: the callers' references at prec + 200
           bits decide all but the exact and near-exact cases) */
        if (arf_equal(lo, hi))
            ok = arf_equal(r, lo);
        else
            ok = (arf_cmp(lo, r) <= 0 && arf_cmp(r, hi) <= 0);

        /* the inexact flag: 0 exactly when the value is exact and fits
           (an exact result the reference does not know to be exact must
           at least lie in the ball) */
        if (ok)
        {
            if (exact_got && !exact_expected)
                ok = arb_contains_arf(b, r);
            else
                ok = (exact_got == (exact_expected && arf_bits(arb_midref(b)) <= prec));
        }
    }
    else
    {
        /* |r - mid| <= rad + tol ulp(r) */
        arf_set_mag(u, arb_radref(b));
        arf_sub(s, r, arb_midref(b), prec + 200, ARF_RND_DOWN);
        arf_abs(s, s);
        if (arf_is_zero(r))
            ok = (arf_cmp(s, u) <= 0);
        else
        {
            arf_t ulp;
            arf_init(ulp);
            arf_set_ui_2exp_si(ulp, 1, -prec);
            arf_mul_2exp_fmpz(ulp, ulp, ARF_EXPREF(r));
            if (rnd == ARF_RND_ACCURATE)
                arf_mul_ui(ulp, ulp, 33, ARF_PREC_EXACT, ARF_RND_DOWN),
                arf_mul_2exp_si(ulp, ulp, -6);    /* 33/64: 0.51 ulp + slack */
            arf_add(u, u, ulp, 64, ARF_RND_UP);
            ok = (arf_cmp(s, u) <= 0);
            arf_clear(ulp);
        }
    }

    if (!ok)
    {
        flint_printf("FAIL: %s, prec = %wd, rnd = %d, inexact = %d\n", name, prec, rnd, inexact);
        flint_printf("x = "); arf_debug(x); flint_printf("\n");
        if (y != NULL) { flint_printf("y = "); arf_debug(y); flint_printf("\n"); }
        flint_printf("r = "); arf_debug(r); flint_printf("\n");
        flint_printf("b = "); arb_printd(b, 40); flint_printf("\n");
        flint_printf("lo = "); arf_debug(lo); flint_printf("\n");
        flint_printf("hi = "); arf_debug(hi); flint_printf("\n");
        flint_abort();
    }

cleanup:
    arf_clear(s);
    arf_clear(lo);
    arf_clear(hi);
    arf_clear(u);
}

static void
arf_elem_verify(const char * name, const arf_t x, const arf_t r, int inexact,
    const arb_t b, slong prec, arf_rnd_t rnd, slong maglim)
{
    arf_elem_verify2(name, x, NULL, r, inexact, b, prec, rnd, maglim);
}

/* Checks r = f(x) (also with aliasing) against the ball g(x) */
static void
arf_elem_check(const char * name, arf_elem_func_t f, arb_elem_func_t g,
    const arf_t x, slong prec, arf_rnd_t rnd, slong maglim, flint_rand_t state)
{
    arf_t r, s;
    arb_t b, bx;
    int inexact, inexact2;

    arf_init(r);
    arf_init(s);
    arb_init(b);
    arb_init(bx);

    inexact = f(r, x, prec, rnd);

    /* aliasing */
    arf_set(s, x);
    inexact2 = f(s, s, prec, rnd);
    if (!arf_equal(r, s) || (inexact != inexact2 && !arf_is_nan(r)))
    {
        flint_printf("FAIL (aliasing): %s, prec = %wd, rnd = %d\n", name, prec, rnd);
        flint_printf("x = "); arf_debug(x); flint_printf("\n");
        flint_printf("r = "); arf_debug(r); flint_printf("\n");
        flint_printf("s = "); arf_debug(s); flint_printf("\n");
        flint_abort();
    }

    arb_set_arf(bx, x);
    g(b, bx, prec + 64);

    arf_elem_verify(name, x, r, inexact, b, prec, rnd, maglim);

    arf_clear(r);
    arf_clear(s);
    arb_clear(b);
    arb_clear(bx);
}

/* the standard loop: random arguments, precisions and modes */
static void
arf_elem_test(const char * name, arf_elem_func_t f, arb_elem_func_t g, int trig, flint_rand_t state, slong iters)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        arf_t x;
        slong prec;
        arf_rnd_t rnd;

        arf_init(x);

        prec = 2 + n_randint(state, 300);
        if (n_randint(state, 8) == 0)
            prec = 2 + n_randint(state, 3000);
        if (n_randint(state, 100) == 0)
            prec = 2 + n_randint(state, 12000);

        rnd = (arf_rnd_t) n_randint(state, 7);

        arf_elem_randtest_arg(x, state, prec);
        arf_elem_check(name, f, g, x, prec, rnd,
            trig ? FLINT_MAX(65536, 4 * prec) : FLINT_MAX(128, 2 * prec), state);

        arf_clear(x);
    }
}

/* the same with the reference computed by a method of the gr context of
   arb (for the functions arb lacks: the generic compositions) */
static void
arf_elem_test_gr(const char * name, arf_elem_func_t f, int method, int trig, flint_rand_t state, slong iters)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        arf_t x, r, s;
        arb_t b, bx;
        gr_ctx_t ctx;
        slong prec, maglim;
        arf_rnd_t rnd;
        int inexact, inexact2;

        arf_init(x);
        arf_init(r);
        arf_init(s);
        arb_init(b);
        arb_init(bx);

        prec = 2 + n_randint(state, 300);
        if (n_randint(state, 8) == 0)
            prec = 2 + n_randint(state, 3000);
        if (n_randint(state, 100) == 0)
            prec = 2 + n_randint(state, 12000);

        rnd = (arf_rnd_t) n_randint(state, 7);
        /* (trig = 2: the limit of exp on x log 10) */
        maglim = (trig == 1) ? FLINT_MAX(65536, 4 * prec) : FLINT_MAX(128, 2 * prec) - 2 * (trig == 2);

        arf_elem_randtest_arg(x, state, prec);

        inexact = f(r, x, prec, rnd);

        arf_set(s, x);
        inexact2 = f(s, s, prec, rnd);
        if (!arf_equal(r, s) || (inexact != inexact2 && !arf_is_nan(r)))
        {
            flint_printf("FAIL (aliasing): %s, prec = %wd, rnd = %d\n", name, prec, rnd);
            flint_printf("x = "); arf_debug(x); flint_printf("\n");
            flint_printf("r = "); arf_debug(r); flint_printf("\n");
            flint_printf("s = "); arf_debug(s); flint_printf("\n");
            flint_abort();
        }

        gr_ctx_init_real_arb(ctx, prec + 64);
        arb_set_arf(bx, x);
        if (((gr_method_unary_op *) ctx->methods)[method](b, bx, ctx) != GR_SUCCESS)
            arb_indeterminate(b);
        gr_ctx_clear(ctx);

        arf_elem_verify(name, x, r, inexact, b, prec, rnd, maglim);

        arf_clear(x);
        arf_clear(r);
        arf_clear(s);
        arb_clear(b);
        arb_clear(bx);
    }
}

#endif
