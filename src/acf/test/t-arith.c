/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acf.h"
#include "acb.h"
#include "fmpz.h"

/* acf_inv, acf_div, acf_sqrt and acf_rsqrt against acb at a higher
   precision: correctly rounded parts in the strict modes (whenever the
   reference ball rounds to one number), within 1 resp. 0.51 ulp in the
   relaxed modes; acf_approx_rsqrt within a few ulps of the magnitude. */

static void
arith_randtest_part(arf_t x, flint_rand_t state, slong prec)
{
    slong kind = n_randint(state, 12);

    arf_randtest(x, state, 1 + n_randint(state, prec), 4);
    if (kind == 0)
        arf_zero(x);
    else if (kind == 1)
        arf_mul_2exp_si(x, x, -(slong) n_randint(state, 600));
    else if (kind == 2)
        arf_mul_2exp_si(x, x, (slong) n_randint(state, 600));
    else if (kind == 3)
        arf_set_si(x, (slong) n_randint(state, 9) - 4);
    else if (kind == 4)
        arf_set_si_2exp_si(x, n_randint(state, 2) ? 1 : -1, (slong) n_randint(state, 20) - 10);
}

static int
arith_ball_rounds_to(arf_t res, const arb_t r, slong prec, arf_rnd_t rnd)
{
    arf_t lo, hi;
    int ok;
    /* (the endpoints rounded outward at the precision of the midpoint and
       more: exact for an exact ball, and no huge exact sum for a radius
       far below the midpoint) */
    slong wp = FLINT_MAX(prec, (slong) arf_bits(arb_midref(r))) + 64;

    arf_init(lo);
    arf_init(hi);
    arb_get_lbound_arf(lo, r, wp);
    arb_get_ubound_arf(hi, r, wp);
    arf_set_round(lo, lo, prec, rnd);
    arf_set_round(hi, hi, prec, rnd);
    ok = arf_equal(lo, hi);
    if (ok)
        arf_set(res, lo);
    arf_clear(lo);
    arf_clear(hi);
    return ok;
}

static int
arith_within_ulps(const arf_t x, const arb_t ref, const arb_t other, slong prec, double bound)
{
    arb_t t;
    mag_t u, v;
    int ok;

    if (arf_is_special(x) && (other == NULL || !arf_is_zero(x)))
        return arf_is_zero(x) && arb_contains_zero(ref);

    arb_init(t);
    mag_init(u);
    mag_init(v);
    arb_set_arf(t, x);
    arb_sub(t, t, ref, prec + 100);
    if (other == NULL)
        arf_get_mag(u, x);
    else
    {
        /* (the relaxed modes: an ulp of the larger part) */
        mag_t w;
        mag_init(w);
        arb_get_mag(u, ref);
        arb_get_mag(w, other);
        mag_max(u, u, w);
        mag_clear(w);
    }
    mag_mul_2exp_si(u, u, -prec + 1);
    mag_mul_ui(u, u, (ulong) (bound * 1024 + 1));
    mag_mul_2exp_si(u, u, -10);
    arf_get_mag_lower(v, arb_midref(t));
    mag_add(u, u, arb_radref(t));
    ok = mag_cmp(v, u) <= 0;
    arb_clear(t);
    mag_clear(u);
    mag_clear(v);
    return ok;
}

static void
arith_check_part(const char * name, const char * part, const arf_t got, const arb_t ref, const arb_t other, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    int ok;

    if (!arb_is_finite(ref))
        return;

    if (rnd < ARF_RND_FAST)
    {
        arf_t r;
        arf_init(r);
        if (arith_ball_rounds_to(r, ref, prec, rnd))
            ok = arf_equal(got, r);
        else
            ok = arith_within_ulps(got, ref, NULL, prec, 1.0);
        arf_clear(r);
    }
    else
        ok = arith_within_ulps(got, ref, other, prec, (rnd == ARF_RND_FAST) ? 1.0 : 0.51);

    if (!ok)
    {
        flint_printf("FAIL: %s, %s part, prec = %wd, rnd = %d\n", name, part, prec, rnd);
        flint_printf("x = "); arf_print(acf_realref(x)); flint_printf(" "); arf_print(acf_imagref(x)); flint_printf("\n");
        if (y != NULL) { flint_printf("y = "); arf_print(acf_realref(y)); flint_printf(" "); arf_print(acf_imagref(y)); flint_printf("\n"); }
        flint_printf("got = "); arf_printd(got, 40); flint_printf("\n");
        flint_printf("ref = "); arb_printd(ref, 40); flint_printf("\n");
        flint_abort();
    }
}

TEST_FUNCTION_START(acf_arith, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        acf_t x, y, r;
        acb_t z, w, zy;
        slong prec, wp, which;
        arf_rnd_t rnd;
        int aliasing;
        const char * name;

        acf_init(x);
        acf_init(y);
        acf_init(r);
        acb_init(z);
        acb_init(w);
        acb_init(zy);

        /* (occasionally beyond 40 limbs: the heap buffers of the kernels) */
        prec = 2 + n_randint(state, (iter % 100 == 0) ? 4000 : (iter % 10 == 0) ? 500 : 140);
        rnd = (arf_rnd_t) n_randint(state, 7);
        which = n_randint(state, 5);
        aliasing = n_randint(state, 4) == 0;

        arith_randtest_part(acf_realref(x), state, prec);
        arith_randtest_part(acf_imagref(x), state, prec);
        arith_randtest_part(acf_realref(y), state, prec);
        arith_randtest_part(acf_imagref(y), state, prec);
        if (n_randint(state, 8) == 0)
        {
            /* a perfect square: x = y^2 */
            acf_mul(x, y, y, ARF_PREC_EXACT, ARF_RND_DOWN);
        }

        arb_set_arf(acb_realref(z), acf_realref(x));
        arb_set_arf(acb_imagref(z), acf_imagref(x));
        arb_set_arf(acb_realref(zy), acf_realref(y));
        arb_set_arf(acb_imagref(zy), acf_imagref(y));
        wp = prec + 64 + 2 * FLINT_BIT_COUNT(prec);

        switch (which)
        {
            case 0:
                name = "inv";
                acb_inv(w, z, wp);
                if (aliasing) { acf_set(r, x); acf_inv(r, r, prec, rnd); } else acf_inv(r, x, prec, rnd);
                break;
            case 1:
                name = "div";
                acb_div(w, z, zy, wp);
                if (aliasing) { acf_set(r, y); acf_div(r, x, r, prec, rnd); } else acf_div(r, x, y, prec, rnd);
                break;
            case 2:
                name = "sqrt";
                acb_sqrt(w, z, wp);
                if (aliasing) { acf_set(r, x); acf_sqrt(r, r, prec, rnd); } else acf_sqrt(r, x, prec, rnd);
                break;
            case 3:
                name = "rsqrt";
                acb_rsqrt(w, z, wp);
                if (aliasing) { acf_set(r, x); acf_rsqrt(r, r, prec, rnd); } else acf_rsqrt(r, x, prec, rnd);
                break;
            default:
                name = "approx_rsqrt";
                acb_rsqrt(w, z, wp);
                acf_approx_rsqrt(r, x, prec, rnd);
                break;
        }

        if (which < 4)
        {
            if (acb_rel_accuracy_bits(w) > prec + 10)
            {
                arith_check_part(name, "real", acf_realref(r), acb_realref(w), acb_imagref(w), x, (which == 1) ? y : NULL, prec, rnd);
                arith_check_part(name, "imag", acf_imagref(r), acb_imagref(w), acb_realref(w), x, (which == 1) ? y : NULL, prec, rnd);
            }
        }
        else if (acb_is_finite(w) && !acb_contains_zero(w))
        {
            /* within 2^(4 - prec) of the magnitude */
            acb_t t;
            mag_t m, u;
            acb_init(t);
            mag_init(m);
            mag_init(u);
            arb_set_arf(acb_realref(t), acf_realref(r));
            arb_set_arf(acb_imagref(t), acf_imagref(r));
            acb_sub(t, t, w, wp);
            acb_get_mag(m, t);
            acb_get_mag(u, w);
            mag_mul_2exp_si(u, u, 4 - prec);
            if (mag_cmp(m, u) > 0)
            {
                flint_printf("FAIL: approx_rsqrt, prec = %wd\n", prec);
                acf_printd(x, 30); flint_printf("\n");
                acf_printd(r, 30); flint_printf("\n");
                acb_printd(w, 30); flint_printf("\n");
                flint_abort();
            }
            acb_clear(t);
            mag_clear(m);
            mag_clear(u);
        }

        acf_clear(x);
        acf_clear(y);
        acf_clear(r);
        acb_clear(z);
        acb_clear(w);
        acb_clear(zy);
    }

    /* parts near a rounding boundary: x = 2^e (m + c 2^-k i) (and y
       likewise) with k far beyond the precision, whose quotients and
       roots lie within about 2^-k or 2^-2k of a dyadic number, against
       acb at a precision resolving that distance; and parts of exponents
       beyond the kernels (2^57 and bignums), against acb */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        acf_t x, y, r;
        acb_t z, zy, w;
        slong prec, k, wp, which;
        arf_rnd_t rnd;
        int i;

        acf_init(x);
        acf_init(y);
        acf_init(r);
        acb_init(z);
        acb_init(zy);
        acb_init(w);

        prec = 2 + n_randint(state, 140);
        k = prec + 8 + n_randint(state, 3000);
        which = n_randint(state, 4);

        for (i = 0; i < 2; i++)
        {
            acf_ptr v = i ? y : x;
            slong m = (n_randint(state, 2) ? 1 : 2 * (slong) n_randint(state, 5) + 1) * (n_randint(state, 2) ? 1 : -1);
            arf_set_si(acf_realref(v), m);
            arf_set_si_2exp_si(acf_imagref(v), (slong) n_randint(state, 9) - 4, -(slong) (n_randint(state, 2) ? k : k / 2 + 1));
            if (n_randint(state, 2))
                arf_swap(acf_realref(v), acf_imagref(v));
            if (n_randint(state, 4) == 0)
            {
                /* huge exponents, the same for both parts */
                fmpz_t e;
                fmpz_init(e);
                fmpz_set_ui(e, n_randint(state, 2) ? (UWORD(1) << 57) : 1);
                if (fmpz_is_one(e))
                    fmpz_mul_2exp(e, e, 70);
                if (n_randint(state, 2))
                    fmpz_neg(e, e);
                arf_mul_2exp_fmpz(acf_realref(v), acf_realref(v), e);
                arf_mul_2exp_fmpz(acf_imagref(v), acf_imagref(v), e);
                fmpz_clear(e);
            }
            else
            {
                slong e = (slong) n_randint(state, 41) - 20;
                arf_mul_2exp_si(acf_realref(v), acf_realref(v), e);
                arf_mul_2exp_si(acf_imagref(v), acf_imagref(v), e);
            }
        }
        if (which >= 2 && n_randint(state, 4) == 0)
            acf_mul(x, y, y, ARF_PREC_EXACT, ARF_RND_DOWN);

        arb_set_arf(acb_realref(z), acf_realref(x));
        arb_set_arf(acb_imagref(z), acf_imagref(x));
        arb_set_arf(acb_realref(zy), acf_realref(y));
        arb_set_arf(acb_imagref(zy), acf_imagref(y));
        wp = 3 * k + 2 * prec + 200;

        if (which == 0) acb_inv(w, z, wp);
        else if (which == 1) acb_div(w, z, zy, wp);
        else if (which == 2) acb_sqrt(w, z, wp);
        else acb_rsqrt(w, z, wp);

        for (rnd = 0; rnd <= ARF_RND_ACCURATE; rnd++)
        {
            const char * name = (which == 0) ? "inv" : (which == 1) ? "div" : (which == 2) ? "sqrt" : "rsqrt";
            if (which == 0) acf_inv(r, x, prec, rnd);
            else if (which == 1) acf_div(r, x, y, prec, rnd);
            else if (which == 2) acf_sqrt(r, x, prec, rnd);
            else acf_rsqrt(r, x, prec, rnd);
            if (acb_rel_accuracy_bits(w) > prec + 10)
            {
                arith_check_part(name, "real (near)", acf_realref(r), acb_realref(w), acb_imagref(w), x, (which == 1) ? y : NULL, prec, rnd);
                arith_check_part(name, "imag (near)", acf_imagref(r), acb_imagref(w), acb_realref(w), x, (which == 1) ? y : NULL, prec, rnd);
            }
        }

        acf_clear(x);
        acf_clear(y);
        acf_clear(r);
        acb_clear(z);
        acb_clear(zy);
        acb_clear(w);
    }

    /* nan and infinite parts: nan results (with a nonzero other part) */
    for (iter = 0; iter < 100; iter++)
    {
        acf_t x, y, r;
        int which = iter % 5, i;

        acf_init(x);
        acf_init(y);
        acf_init(r);
        for (i = 0; i < 2; i++)
        {
            arf_ptr v = i ? acf_imagref(x) : acf_realref(x);
            switch (n_randint(state, 3))
            {
                case 0: arf_nan(v); break;
                case 1: arf_pos_inf(v); break;
                default: arf_neg_inf(v); break;
            }
        }
        arf_set_si(n_randint(state, 2) ? acf_realref(x) : acf_imagref(x), (slong) n_randint(state, 7) - 3);
        if (arf_is_zero(acf_realref(x)) || arf_is_zero(acf_imagref(x)))
            arf_set_si(arf_is_zero(acf_realref(x)) ? acf_realref(x) : acf_imagref(x), 5);
        arf_set_si(acf_realref(y), 3);
        arf_set_si(acf_imagref(y), -1);
        if (which == 0) acf_inv(r, x, 64, ARF_RND_DOWN);
        else if (which == 1) acf_div(r, x, y, 64, ARF_RND_DOWN);
        else if (which == 2) acf_div(r, y, x, 64, ARF_RND_DOWN);
        else if (which == 3) acf_sqrt(r, x, 64, ARF_RND_DOWN);
        else acf_rsqrt(r, x, 64, ARF_RND_DOWN);
        if (!arf_is_nan(acf_realref(r)) || !arf_is_nan(acf_imagref(r)))
        {
            flint_printf("FAIL: nonfinite (%d)\n", which);
            acf_printd(x, 10); flint_printf("\n");
            acf_printd(r, 10); flint_printf("\n");
            flint_abort();
        }
        acf_clear(x);
        acf_clear(y);
        acf_clear(r);
    }

    /* far beyond any working precision: 1 + eps i with eps = 2^-(2^62)
       (a bignum exponent), whose reciprocal has the real part 1 / (1 +
       eps^2), just below 1, and the imaginary part -eps / (1 + eps^2),
       just above -eps; the square root (1 + eps^2/8 + ..., eps/2 (1 -
       ...)), each part just beyond a power of two */
    {
        acf_t x, r;
        arf_t eps, e1, e2, t;
        fmpz_t e;
        int rnd;

        acf_init(x);
        acf_init(r);
        arf_init(eps);
        arf_init(e1);
        arf_init(e2);
        arf_init(t);
        fmpz_init(e);

        fmpz_one(e);
        fmpz_mul_2exp(e, e, 62);
        fmpz_neg(e, e);
        arf_one(eps);
        arf_mul_2exp_fmpz(eps, eps, e);
        arf_one(acf_realref(x));
        arf_set(acf_imagref(x), eps);

        for (rnd = 0; rnd <= ARF_RND_NEAR; rnd++)
        {
            int st, toward0 = (rnd == ARF_RND_DOWN), up = (rnd == ARF_RND_UP);

            /* inv: the real part just below 1, the imaginary part just
               above -eps (the numbers beneath a power of two 2^t at 64
               bits: 2^t (1 - 2^-64)) */
            st = acf_inv(r, x, 64, rnd);
            arf_one(e1);
            if (toward0 || rnd == ARF_RND_FLOOR)
            {
                arf_set_si_2exp_si(t, 1, -64);
                arf_sub(e1, e1, t, ARF_PREC_EXACT, ARF_RND_DOWN);
            }
            arf_neg(e2, eps);
            if (toward0 || rnd == ARF_RND_CEIL)
            {
                arf_mul_2exp_si(t, e2, -64);
                arf_sub(e2, e2, t, ARF_PREC_EXACT, ARF_RND_DOWN);
            }
            if (!arf_equal(acf_realref(r), e1) || !arf_equal(acf_imagref(r), e2) || st != 3)
            {
                flint_printf("FAIL: inv(1 + eps i), rnd = %d\n", rnd);
                acf_printd(r, 20); flint_printf("\n");
                flint_abort();
            }

            /* sqrt: the real part just above 1, the imaginary part just
               below eps/2 */
            st = acf_sqrt(r, x, 64, rnd);
            arf_one(e1);
            if (up || rnd == ARF_RND_CEIL)
            {
                arf_set_si_2exp_si(t, 1, -63);
                arf_add(e1, e1, t, ARF_PREC_EXACT, ARF_RND_DOWN);
            }
            arf_mul_2exp_si(e2, eps, -1);
            if (toward0 || rnd == ARF_RND_FLOOR)
            {
                arf_mul_2exp_si(t, e2, -64);
                arf_sub(e2, e2, t, ARF_PREC_EXACT, ARF_RND_DOWN);
            }
            if (!arf_equal(acf_realref(r), e1) || !arf_equal(acf_imagref(r), e2) || st != 3)
            {
                flint_printf("FAIL: sqrt(1 + eps i), rnd = %d\n", rnd);
                acf_printd(r, 20); flint_printf("\n");
                flint_abort();
            }
        }

        acf_clear(x);
        acf_clear(r);
        arf_clear(eps);
        arf_clear(e1);
        arf_clear(e2);
        arf_clear(t);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
