/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "arb.h"
#include "mp_real.h"

/* mp_real_log1p_bits, expm1, sinh_cosh, tanh, asin, acos, asinh, acosh,
   atanh: the output ball contains f at a random point of the input ball,
   for tiny, unit and large arguments, arguments near +-1, small integers
   and powers of two, inexact arguments and aliased outputs; exact
   arguments give a relative accuracy of about prec bits; the domain
   reports are consistent with the input ball (a failure is allowed only
   when the ball comes within 2^31 radii of a boundary of the domain). */

#define F_LOG1P 0
#define F_EXPM1 1
#define F_SINH 2
#define F_COSH 3
#define F_TANH 4
#define F_ASIN 5
#define F_ACOS 6
#define F_ASINH 7
#define F_ACOSH 8
#define F_ATANH 9
#define NUM_F 10

static const char * fname[] = { "log1p", "expm1", "sinh", "cosh", "tanh",
    "asin", "acos", "asinh", "acosh", "atanh" };

static int
eval(mp_real_t y, const mp_real_t x, int f, slong prec)
{
    switch (f)
    {
        case F_LOG1P: return mp_real_log1p_bits(y, x, prec);
        case F_EXPM1: mp_real_expm1_bits(y, x, prec); return 1;
        case F_SINH: mp_real_sinh_cosh_bits(y, NULL, x, prec); return 1;
        case F_COSH: mp_real_sinh_cosh_bits(NULL, y, x, prec); return 1;
        case F_TANH: mp_real_tanh_bits(y, x, prec); return 1;
        case F_ASIN: return mp_real_asin_bits(y, x, prec);
        case F_ACOS: return mp_real_acos_bits(y, x, prec);
        case F_ASINH: mp_real_asinh_bits(y, x, prec); return 1;
        case F_ACOSH: return mp_real_acosh_bits(y, x, prec);
        default: return mp_real_atanh_bits(y, x, prec);
    }
}

static void
ref(arb_t y, const arb_t x, int f, slong prec)
{
    switch (f)
    {
        case F_LOG1P: arb_log1p(y, x, prec); break;
        case F_EXPM1: arb_expm1(y, x, prec); break;
        case F_SINH: arb_sinh(y, x, prec); break;
        case F_COSH: arb_cosh(y, x, prec); break;
        case F_TANH: arb_tanh(y, x, prec); break;
        case F_ASIN: arb_asin(y, x, prec); break;
        case F_ACOS: arb_acos(y, x, prec); break;
        case F_ASINH: arb_asinh(y, x, prec); break;
        case F_ACOSH: arb_acosh(y, x, prec); break;
        default: arb_atanh(y, x, prec); break;
    }
}

/* whether a failure is acceptable: the ball [m +- r] reaches the
   domain boundary b within 2^31 r (or crosses it) */
static int
near_boundary(const arb_t x, int f)
{
    arb_t d, t;
    mag_t r;
    int res;

    arb_init(d);
    arb_init(t);
    mag_init(r);

    /* d = the distance of the midpoint to the boundary on the side of
       the domain (negative outside) */
    if (f == F_LOG1P)
    {
        arb_set_arf(d, arb_midref(x));
        arb_add_ui(d, d, 1, 1000);
    }
    else if (f == F_ACOSH)
    {
        arb_set_arf(d, arb_midref(x));
        arb_sub_ui(d, d, 1, 1000);
    }
    else
    {
        arb_set_arf(d, arb_midref(x));
        arb_abs(d, d);
        arb_sub_ui(d, d, 1, 1000);
        arb_neg(d, d);
    }

    mag_mul_2exp_si(r, arb_radref(x), 31);
    arb_zero(t);
    arf_set_mag(arb_midref(t), r);
    res = !arb_gt(d, t);

    arb_clear(d);
    arb_clear(t);
    mag_clear(r);
    return res;
}

TEST_FUNCTION_START(mp_real_elem_composed, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        int f = n_randint(state, NUM_F), kind, inexact, alias, ok, expected;
        slong prec, L, e2, wp;
        mp_real_t x, y, d;
        arb_t xa, ya, t, ty;
        arf_t pt;
        nn_ptr p;

        if (iter % 100 == 0)
            prec = 2 + n_randint(state, 5000);
        else if (iter % 10 == 0)
            prec = 2 + n_randint(state, 1000);
        else
            prec = 2 + n_randint(state, 300);

        kind = n_randint(state, 8);

        mp_real_init(x);
        mp_real_init(y);
        mp_real_init(d);
        arb_init(xa);
        arb_init(ya);
        arb_init(t);
        arb_init(ty);
        arf_init(pt);

        L = 1 + n_randint(state, prec / FLINT_BITS + 3);
        p = flint_malloc(L * sizeof(ulong));
        flint_mpn_rrandom(p, state, L);
        if (p[L - 1] == 0)
            p[L - 1] = 1;

        switch (kind)
        {
            case 0: e2 = -(slong) n_randint(state, 2 * prec + 200); break;
            case 1: e2 = -(slong) n_randint(state, 64); break;
            case 2: e2 = 1 + n_randint(state, 20); break;
            case 3: e2 = (f == F_EXPM1 || f == F_SINH || f == F_COSH)
                        ? (slong) n_randint(state, FLINT_BITS - 6)   /* (exp's range) */
                        : (slong) n_randint(state, 3000) - 1500; break;
            default: e2 = (slong) n_randint(state, 8) - 4; break;
        }
        _mp_real_set_mpn_2exp(x, p, L, e2 - FLINT_BITS * L);

        if (kind == 4 || kind == 5)
        {
            /* 1 +- 2^-k (+ a perturbation below) */
            mp_real_set_ui(d, 1);
            mp_real_mul_2exp_si(d, d, -(slong) n_randint(state,
                (kind == 4) ? prec + 100 : 64) - 1);
            mp_real_set(y, x);
            mp_real_set_ui(x, 1);
            if (n_randint(state, 2))
                mp_real_add(x, x, d, 10000);
            else
                mp_real_sub(x, x, d, 10000);
            if (n_randint(state, 2))
            {
                mp_real_mul_2exp_si(y, y, -(slong) n_randint(state, prec + 300) - 2);
                mp_real_add(x, x, y, 10000);
            }
            x->err = 0;
        }

        if (kind == 6)
        {
            mp_real_set_ui(x, n_randint(state, 5));
            mp_real_mul_2exp_si(x, x, (slong) n_randint(state, 6) - 3);
        }

        if (n_randint(state, 2) && !(f == F_ACOSH && n_randint(state, 8) != 0))
            mp_real_neg(x, x);
        if (f == F_LOG1P && kind == 4 && n_randint(state, 2))
        {
            /* near -1 */
            mp_real_set_ui(d, 2);
            mp_real_sub(x, x, d, 10000);
        }

        inexact = (n_randint(state, 4) == 0);
        if (inexact)
        {
            if (mp_real_is_zero(x))
                mp_real_add_error_2exp_si(x, -(slong) n_randint(state, 100));
            else
                mp_real_add_error_2exp_si(x, mp_real_abs_bound_lt_2exp_si(x)
                    - 1 - (slong) n_randint(state, prec + 20));
        }
        mp_real_get_arb(xa, x);

        alias = (n_randint(state, 4) == 0);
        if (alias)
        {
            mp_real_set(y, x);
            ok = eval(y, y, f, prec);
        }
        else
            ok = eval(y, x, f, prec);
        mp_real_get_arb(ya, y);

        /* the domain: inside (for exact arguments, the closed domain
           of asin, acos and acosh) must succeed unless near the
           boundary, outside must fail; with the exact endpoints lo, hi
           of the input ball (arb's radius is rounded up) */
        {
            int inside, near = 0;
            arb_t lo, hi, rr;
            mp_real_t xm, xr;

            arb_init(lo);
            arb_init(hi);
            arb_init(rr);
            mp_real_init(xm);
            mp_real_init(xr);
            mp_real_set(xm, x);
            xm->err = 0;
            mp_real_get_arb(lo, xm);
            if (x->err != 0)
            {
                ulong xe = x->err;
                _mp_real_set_mpn_2exp(xr, &xe, 1, FLINT_BITS * (x->exp - x->size));
                mp_real_get_arb(rr, xr);
            }
            arb_add(hi, lo, rr, ARF_PREC_EXACT);
            arb_sub(lo, lo, rr, ARF_PREC_EXACT);

            switch (f)
            {
                case F_LOG1P:
                    inside = (arf_cmp_si(arb_midref(lo), -1) > 0);
                    near = !arb_is_exact(xa) && near_boundary(xa, f);
                    break;
                case F_ASIN: case F_ACOS:
                    if (x->err == 0)
                        inside = (arf_cmpabs_ui(arb_midref(lo), 1) <= 0);
                    else
                        inside = (arf_cmp_si(arb_midref(lo), -1) > 0 && arf_cmp_si(arb_midref(hi), 1) < 0);
                    near = !arb_is_exact(xa) && near_boundary(xa, f);
                    break;
                case F_ACOSH:
                    if (x->err == 0)
                        inside = (arf_cmp_si(arb_midref(lo), 1) >= 0);
                    else
                        inside = (arf_cmp_si(arb_midref(lo), 1) > 0);
                    near = !arb_is_exact(xa) && near_boundary(xa, f);
                    break;
                case F_ATANH:
                    inside = (arf_cmp_si(arb_midref(lo), -1) > 0 && arf_cmp_si(arb_midref(hi), 1) < 0);
                    near = !arb_is_exact(xa) && near_boundary(xa, f);
                    break;
                default:
                    inside = 1;
            }
            (void) expected;

            arb_clear(lo);
            arb_clear(hi);
            arb_clear(rr);
            mp_real_clear(xm);
            mp_real_clear(xr);

            if ((!inside && ok) || (inside && !near && !ok))
            {
                flint_printf("FAIL: domain, %s, ok = %d\nx = ", fname[f], ok);
                arb_printd(xa, 30);
                flint_printf("\n");
                flint_abort();
            }
        }

        if (ok)
        {
            /* a random point of the input ball */
            arf_set(pt, arb_midref(xa));
            if (inexact)
            {
                arf_t r;
                arf_init(r);
                arf_set_mag(r, arb_radref(xa));
                arf_mul_2exp_si(r, r, -1);
                if (n_randint(state, 2))
                    arf_neg(r, r);
                arf_add(pt, pt, r, ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_clear(r);
            }
            arb_set_arf(t, pt);
            wp = 2 * prec + 100 + (arf_is_zero(pt) ? 0
                : FLINT_ABS(fmpz_get_si(ARF_EXPREF(pt))));
            ref(ty, t, f, wp);

            if (!arb_overlaps(ty, ya))
            {
                flint_printf("FAIL: containment, %s, kind %d, prec %wd, "
                    "inexact %d, alias %d\nx = ", fname[f], kind, prec, inexact, alias);
                arb_printd(xa, 30);
                flint_printf("\ny = "); arb_printd(ya, 30);
                flint_printf("\n    "); arb_printd(ty, 30);
                flint_printf("\n");
                flint_abort();
            }

            /* for inexact arguments, the radius propagated about as well
               as by arb */
            if (inexact && !arb_is_exact(ya))
            {
                arb_t ra;
                mag_t bound, u;

                arb_init(ra);
                mag_init(bound);
                mag_init(u);
                ref(ra, xa, f, wp);

                if (arb_is_finite(ra))
                {
                    mag_mul_2exp_si(bound, arb_radref(ra), 10);
                    arb_get_mag(u, ra);
                    mag_mul_2exp_si(u, u, -prec + 10);
                    mag_add(bound, bound, u);

                    if (mag_cmp(arb_radref(ya), bound) > 0)
                    {
                        flint_printf("FAIL: radius, %s, kind %d, prec %wd\nx = ", fname[f], kind, prec);
                        arb_printd(xa, 30);
                        flint_printf("\ny = "); arb_printd(ya, 30);
                        flint_printf("\n    "); arb_printd(ra, 30);
                        flint_printf("\n");
                        flint_abort();
                    }
                }

                arb_clear(ra);
                mag_clear(bound);
                mag_clear(u);
            }

            if (!inexact && !arb_is_exact(ya)
                && arb_rel_accuracy_bits(ya) < prec - 4)
            {
                flint_printf("FAIL: accuracy, %s, kind %d, prec %wd, acc %wd\nx = ",
                    fname[f], kind, prec, arb_rel_accuracy_bits(ya));
                arb_printd(xa, 30);
                flint_printf("\ny = "); arb_printd(ya, 30);
                flint_printf("\n");
                flint_abort();
            }
        }

        flint_free(p);
        mp_real_clear(x);
        mp_real_clear(y);
        mp_real_clear(d);
        arb_clear(xa);
        arb_clear(ya);
        arb_clear(t);
        arb_clear(ty);
        arf_clear(pt);
    }

    /* tiny arguments, down to a quarter of the safe exponent range: the
       remainder bounds must not pad the outputs down to their scale
       (which would take up to B^(2^55) limbs; t-tiny_arg.c tests the
       functions these build on) */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        mp_real_t x, y, z;
        slong prec = 2 + n_randint(state, 500), e, maxsize;
        ulong d = n_randtest_not_zero(state);
        int f;

        mp_real_init(x);
        mp_real_init(y);
        mp_real_init(z);

        e = -(slong) n_randint(state, FLINT_BITS * (MP_REAL_EXP_MAX / 4));
        _mp_real_set_mpn_2exp(x, &d, 1, e);
        if (n_randint(state, 2))
            mp_real_neg(x, x);
        if (n_randint(state, 2))
            mp_real_add_error_2exp_si(x, mp_real_abs_bound_lt_2exp_si(x) - prec - n_randint(state, 100));

        maxsize = mp_real_prec_bits(prec) + 8;

        for (f = 0; f < NUM_F; f++)
        {
            mp_real_zero(z);
            if (f == F_SINH)
                mp_real_sinh_cosh_bits(y, z, x, prec);
            else if (f != F_COSH && f != F_ACOSH)
                eval(y, x, f, prec);

            if (y->size > maxsize || z->size > maxsize)
            {
                flint_printf("FAIL: tiny argument, %s, e = %wd, prec = %wd, sizes %wd %wd\n",
                    fname[f], e, prec, y->size, z->size);
                flint_abort();
            }
        }

        mp_real_clear(x);
        mp_real_clear(y);
        mp_real_clear(z);
    }

    TEST_FUNCTION_END(state);
}
