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

/* The elementary functions of acf against acb at a higher precision:
   in the strict modes each part is the correctly rounded value (checked
   whenever the reference ball rounds to one number), in ARF_RND_FAST
   within an ulp and in ARF_RND_ACCURATE within 0.51 ulp of the value. */

typedef int (*acf_func_t)(acf_t, const acf_t, slong, arf_rnd_t);
typedef void (*acb_func_t)(acb_t, const acb_t, slong);

static int acf_exp_w(acf_t r, const acf_t x, slong p, arf_rnd_t rnd) { return acf_exp(r, x, p, rnd); }
static int acf_sinh_w(acf_t r, const acf_t x, slong p, arf_rnd_t rnd) { return acf_sinh(r, x, p, rnd); }
static int acf_cosh_w(acf_t r, const acf_t x, slong p, arf_rnd_t rnd) { return acf_cosh(r, x, p, rnd); }

/* the references missing from acb: the reciprocal functions as f(1/z) and
   the pi versions of the inverse functions as f(z) / pi */
static void acb_sec_pi_w(acb_t r, const acb_t z, slong p) { acb_cos_pi(r, z, p); acb_inv(r, r, p); }
static void acb_acot_w(acb_t r, const acb_t z, slong p) { acb_inv(r, z, p); acb_atan(r, r, p); }
static void acb_asec_w(acb_t r, const acb_t z, slong p) { acb_inv(r, z, p); acb_acos(r, r, p); }
static void acb_acsc_w(acb_t r, const acb_t z, slong p) { acb_inv(r, z, p); acb_asin(r, r, p); }
static void acb_acoth_w(acb_t r, const acb_t z, slong p) { acb_inv(r, z, p); acb_atanh(r, r, p); }
static void acb_asech_w(acb_t r, const acb_t z, slong p) { acb_inv(r, z, p); acb_acosh(r, r, p); }
static void acb_acsch_w(acb_t r, const acb_t z, slong p) { acb_inv(r, z, p); acb_asinh(r, r, p); }
static void acb_div_pi(acb_t r, slong p) { arb_t t; arb_init(t); arb_const_pi(t, p); acb_div_arb(r, r, t, p); arb_clear(t); }
static void acb_asin_pi_w(acb_t r, const acb_t z, slong p) { acb_asin(r, z, p); acb_div_pi(r, p); }
static void acb_acos_pi_w(acb_t r, const acb_t z, slong p) { acb_acos(r, z, p); acb_div_pi(r, p); }
static void acb_atan_pi_w(acb_t r, const acb_t z, slong p) { acb_atan(r, z, p); acb_div_pi(r, p); }
static void acb_tan_w(acb_t r, const acb_t z, slong p) { acb_tan(r, z, p); }
static void acb_cot_w(acb_t r, const acb_t z, slong p) { acb_cot(r, z, p); }
static void acb_sec_w(acb_t r, const acb_t z, slong p) { acb_sec(r, z, p); }
static void acb_csc_w(acb_t r, const acb_t z, slong p) { acb_csc(r, z, p); }
static void acb_tanh_w(acb_t r, const acb_t z, slong p) { acb_tanh(r, z, p); }
static void acb_coth_w(acb_t r, const acb_t z, slong p) { acb_coth(r, z, p); }
static void acb_sech_w(acb_t r, const acb_t z, slong p) { acb_sech(r, z, p); }
static void acb_csch_w(acb_t r, const acb_t z, slong p) { acb_csch(r, z, p); }

static const struct { const char * name; acf_func_t f; acb_func_t g; } funcs[] = {
    { "exp", acf_exp_w, acb_exp }, { "exp_pi_i", acf_exp_pi_i, acb_exp_pi_i },
    { "log", acf_log, acb_log }, { "sin", acf_sin, acb_sin }, { "cos", acf_cos, acb_cos },
    { "sin_pi", acf_sin_pi, acb_sin_pi }, { "cos_pi", acf_cos_pi, acb_cos_pi },
    { "sinh", acf_sinh_w, acb_sinh }, { "cosh", acf_cosh_w, acb_cosh },
    { "tan", acf_tan, acb_tan_w }, { "cot", acf_cot, acb_cot_w }, { "sec", acf_sec, acb_sec_w }, { "csc", acf_csc, acb_csc_w },
    { "tan_pi", acf_tan_pi, acb_tan_pi }, { "cot_pi", acf_cot_pi, acb_cot_pi }, { "sec_pi", acf_sec_pi, acb_sec_pi_w }, { "csc_pi", acf_csc_pi, acb_csc_pi },
    { "tanh", acf_tanh, acb_tanh_w }, { "coth", acf_coth, acb_coth_w }, { "sech", acf_sech, acb_sech_w }, { "csch", acf_csch, acb_csch_w },
    { "expm1", acf_expm1, acb_expm1 }, { "log1p", acf_log1p, acb_log1p }, { "sinc", acf_sinc, acb_sinc }, { "sinc_pi", acf_sinc_pi, acb_sinc_pi },
    { "asin", acf_asin, acb_asin }, { "acos", acf_acos, acb_acos }, { "atan", acf_atan, acb_atan },
    { "acot", acf_acot, acb_acot_w }, { "asec", acf_asec, acb_asec_w }, { "acsc", acf_acsc, acb_acsc_w },
    { "asinh", acf_asinh, acb_asinh }, { "acosh", acf_acosh, acb_acosh }, { "atanh", acf_atanh, acb_atanh },
    { "acoth", acf_acoth, acb_acoth_w }, { "asech", acf_asech, acb_asech_w }, { "acsch", acf_acsch, acb_acsch_w },
    { "asin_pi", acf_asin_pi, acb_asin_pi_w }, { "acos_pi", acf_acos_pi, acb_acos_pi_w }, { "atan_pi", acf_atan_pi, acb_atan_pi_w },
};

#define NF ((slong) (sizeof(funcs) / sizeof(funcs[0])))

static void
randtest_part(arf_t x, flint_rand_t state, slong prec)
{
    slong kind = n_randint(state, 16);

    arf_randtest(x, state, 1 + n_randint(state, prec), 4);
    if (kind == 0)
        arf_zero(x);
    else if (kind == 1)
        arf_mul_2exp_si(x, x, -(slong) n_randint(state, 400));
    else if (kind == 5)
        arf_mul_2exp_si(x, x, -(slong) n_randint(state, 3000));
    else if (kind == 2)
        arf_mul_2exp_si(x, x, (slong) n_randint(state, 60));
    else if (kind == 3)
        arf_set_si(x, (slong) n_randint(state, 9) - 4);
    else if (kind == 4)
    {
        /* near the unit circle for log */
        arf_set_ui(x, 1);
        if (n_randint(state, 2))
        {
            arf_t t;
            arf_init(t);
            arf_randtest(t, state, prec, 0);
            arf_mul_2exp_si(t, t, -(slong) n_randint(state, 2 * prec + 10));
            arf_add(x, x, t, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_clear(t);
        }
    }
}

/* whether the value V of the ball r rounds to one number at prec bits
   in rnd (then written to res) */
static int
ball_rounds_to(arf_t res, const arb_t r, slong prec, arf_rnd_t rnd)
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

/* |x - V| <= bound ulps of x (plus the radius of the reference) */
static int
within_ulps(const arf_t x, const arb_t ref, const arb_t other, slong prec, double bound)
{
    arb_t t;
    mag_t u;
    int ok;

    if (arf_is_special(x) && (other == NULL || !arf_is_zero(x)))
        return arf_is_zero(x) && arb_contains_zero(ref);

    arb_init(t);
    mag_init(u);
    arb_set_arf(t, x);
    arb_sub(t, t, ref, prec + 100);
    /* u = bound ulps of x */
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
    {
        /* |mid(t)| <= u + rad(t) */
        mag_t v;
        mag_init(v);
        arf_get_mag_lower(v, arb_midref(t));
        mag_add(u, u, arb_radref(t));
        ok = mag_cmp(v, u) <= 0;
        mag_clear(v);
    }
    arb_clear(t);
    mag_clear(u);
    return ok;
}

static void
check_part(const char * name, const char * part, const arf_t got, const arb_t ref, const arb_t other, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    int ok;

    if (!arb_is_finite(ref))
        return;

    if (rnd < ARF_RND_FAST)
    {
        arf_t r;
        arf_init(r);
        if (ball_rounds_to(r, ref, prec, rnd))
            ok = arf_equal(got, r);
        else
            ok = within_ulps(got, ref, NULL, prec, 1.0);
        arf_clear(r);
    }
    else
        ok = within_ulps(got, ref, other, prec, (rnd == ARF_RND_FAST) ? 1.0 : 0.51);

    if (!ok)
    {
        flint_printf("FAIL: %s, %s part, prec = %wd, rnd = %d\n", name, part, prec, rnd);
        flint_printf("x = "); acf_printd(x, 30); flint_printf("\n");
        flint_printf("x = "); arf_print(acf_realref(x)); flint_printf(" "); arf_print(acf_imagref(x)); flint_printf("\n");
        if (y != NULL) { flint_printf("y = "); acf_printd(y, 30); flint_printf("\n"); }
        if (y != NULL) { flint_printf("y = "); arf_print(acf_realref(y)); flint_printf(" "); arf_print(acf_imagref(y)); flint_printf("\n"); }
        flint_printf("got = "); arf_printd(got, 40); flint_printf("\n");
        flint_printf("ref = "); arb_printd(ref, 40); flint_printf("\n");
        flint_abort();
    }
}

TEST_FUNCTION_START(acf_elem, state)
{
    slong iter;

    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        acf_t x, y, r, r2, xcopy;
        acb_t z, w, w2;
        slong prec, wp, which;
        arf_rnd_t rnd;
        int aliasing;

        acf_init(x);
        acf_init(y);
        acf_init(r);
        acf_init(r2);
        acf_init(xcopy);
        acb_init(z);
        acb_init(w);
        acb_init(w2);

        prec = 2 + n_randint(state, (iter % 10 == 0) ? 500 : 140);
        rnd = (arf_rnd_t) n_randint(state, 7);
        which = n_randint(state, NF + 4);
        aliasing = n_randint(state, 4) == 0;

        randtest_part(acf_realref(x), state, prec);
        randtest_part(acf_imagref(x), state, prec);
        randtest_part(acf_realref(y), state, prec);
        randtest_part(acf_imagref(y), state, prec);
        if (n_randint(state, 3) == 0)
        {
            arf_mul_2exp_si(acf_realref(y), acf_realref(y), -5 - (slong) n_randint(state, 30));
            arf_mul_2exp_si(acf_imagref(y), acf_imagref(y), -5 - (slong) n_randint(state, 30));
        }
        acf_set(xcopy, x);

        arb_set_arf(acb_realref(z), acf_realref(x));
        arb_set_arf(acb_imagref(z), acf_imagref(x));
        wp = prec + 64 + 2 * FLINT_BIT_COUNT(FLINT_ABS(arf_is_zero(acf_realref(x)) ? 0 : fmpz_get_si(ARF_EXPREF(acf_realref(x)))))
            + 2 * FLINT_BIT_COUNT(FLINT_ABS(arf_is_zero(acf_imagref(x)) ? 0 : fmpz_get_si(ARF_EXPREF(acf_imagref(x)))));

        if (which < NF)
        {
            /* (acb's accuracy for the smaller part of the value drops with
               the gap between the parts of the argument) */
            if (!arf_is_zero(acf_realref(x)) && !arf_is_zero(acf_imagref(x)) && arf_is_finite(acf_realref(x)) && arf_is_finite(acf_imagref(x)))
            {
                slong ga = fmpz_get_si(ARF_EXPREF(acf_realref(x))), gb = fmpz_get_si(ARF_EXPREF(acf_imagref(x)));
                wp += FLINT_ABS(ga - gb);
            }

            funcs[which].g(w, z, wp);
            if (aliasing)
            {
                acf_set(r, x);
                funcs[which].f(r, r, prec, rnd);
            }
            else
                funcs[which].f(r, x, prec, rnd);
            /* (poles: acf gives nan parts, acb may have a finite one) */
            if (acb_is_finite(w))
            {
                check_part(funcs[which].name, "real", acf_realref(r), acb_realref(w), acb_imagref(w), x, NULL, prec, rnd);
                check_part(funcs[which].name, "imag", acf_imagref(r), acb_imagref(w), acb_realref(w), x, NULL, prec, rnd);
            }
        }
        else if (which < NF + 3)
        {
            const char * name = (which == NF) ? "sin_cos" : (which == NF + 1) ? "sin_cos_pi" : "sinh_cosh";
            if (which == NF) acb_sin_cos(w, w2, z, wp);
            else if (which == NF + 1) acb_sin_cos_pi(w, w2, z, wp);
            else acb_sinh_cosh(w, w2, z, wp);
            if (aliasing)
            {
                acf_set(r, x);
                if (which == NF) acf_sin_cos(r, r2, r, prec, rnd);
                else if (which == NF + 1) acf_sin_cos_pi(r, r2, r, prec, rnd);
                else acf_sinh_cosh(r, r2, r, prec, rnd);
            }
            else
            {
                if (which == NF) acf_sin_cos(r, r2, x, prec, rnd);
                else if (which == NF + 1) acf_sin_cos_pi(r, r2, x, prec, rnd);
                else acf_sinh_cosh(r, r2, x, prec, rnd);
            }
            check_part(name, "sin real", acf_realref(r), acb_realref(w), acb_imagref(w), x, NULL, prec, rnd);
            check_part(name, "sin imag", acf_imagref(r), acb_imagref(w), acb_realref(w), x, NULL, prec, rnd);
            check_part(name, "cos real", acf_realref(r2), acb_realref(w2), acb_imagref(w2), x, NULL, prec, rnd);
            check_part(name, "cos imag", acf_imagref(r2), acb_imagref(w2), acb_realref(w2), x, NULL, prec, rnd);
        }
        else
        {
            acb_t zy;
            acb_init(zy);
            arb_set_arf(acb_realref(zy), acf_realref(y));
            arb_set_arf(acb_imagref(zy), acf_imagref(y));
            wp += 2 * FLINT_BIT_COUNT(FLINT_ABS(arf_is_zero(acf_realref(y)) ? 0 : fmpz_get_si(ARF_EXPREF(acf_realref(y)))))
                + 2 * FLINT_BIT_COUNT(FLINT_ABS(arf_is_zero(acf_imagref(y)) ? 0 : fmpz_get_si(ARF_EXPREF(acf_imagref(y)))));
            acb_pow(w, z, zy, wp);
            if (aliasing)
            {
                acf_set(r, x);
                acf_pow(r, r, y, prec, rnd);
            }
            else
                acf_pow(r, x, y, prec, rnd);
            /* (acb_pow may return a wide ball for huge |y log x|) */
            if (acb_rel_accuracy_bits(w) > prec + 10 || !acb_is_finite(w))
            {
                check_part("pow", "real", acf_realref(r), acb_realref(w), acb_imagref(w), x, y, prec, rnd);
                check_part("pow", "imag", acf_imagref(r), acb_imagref(w), acb_realref(w), x, y, prec, rnd);
            }
            acb_clear(zy);
        }

        if (!acf_equal(x, xcopy) && !aliasing)
        {
            flint_printf("FAIL: input modified\n");
            flint_abort();
        }

        acf_clear(x);
        acf_clear(y);
        acf_clear(r);
        acf_clear(r2);
        acf_clear(xcopy);
        acb_clear(z);
        acb_clear(w);
        acb_clear(w2);
    }

    /* the register path of pow with a far smaller real part of x and
       of y (shifts of the terms y log x by several limbs), in every
       mode: x = -291173 2^-223 - 123 2^5 i, y = -57 2^-198 - i, and
       random parts of such sizes */
    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        acf_t x, y, r;
        acb_t z, zy, w;
        arf_rnd_t rnd;

        acf_init(x);
        acf_init(y);
        acf_init(r);
        acb_init(z);
        acb_init(zy);
        acb_init(w);

        if (iter == 0)
        {
            arf_set_si_2exp_si(acf_realref(x), -291173, -223);
            arf_set_si_2exp_si(acf_imagref(x), -123, 5);
            arf_set_si_2exp_si(acf_realref(y), -57, -198);
            arf_set_si_2exp_si(acf_imagref(y), -1, 0);
        }
        else
        {
            arf_set_si_2exp_si(acf_realref(x), (slong) n_randint(state, 1000000) + 1, -(slong) n_randint(state, 400));
            arf_set_si_2exp_si(acf_imagref(x), (slong) n_randint(state, 1000) + 1, (slong) n_randint(state, 20) - 10);
            arf_set_si_2exp_si(acf_realref(y), (slong) n_randint(state, 1000) + 1, -(slong) n_randint(state, 400));
            arf_set_si_2exp_si(acf_imagref(y), (slong) n_randint(state, 7) + 1, (slong) n_randint(state, 4) - 2);
            if (n_randint(state, 2)) arf_neg(acf_realref(x), acf_realref(x));
            if (n_randint(state, 2)) arf_neg(acf_imagref(x), acf_imagref(x));
            if (n_randint(state, 2)) arf_neg(acf_realref(y), acf_realref(y));
            if (n_randint(state, 2)) arf_neg(acf_imagref(y), acf_imagref(y));
            if (n_randint(state, 2)) arf_swap(acf_realref(x), acf_imagref(x));
        }

        arb_set_arf(acb_realref(z), acf_realref(x));
        arb_set_arf(acb_imagref(z), acf_imagref(x));
        arb_set_arf(acb_realref(zy), acf_realref(y));
        arb_set_arf(acb_imagref(zy), acf_imagref(y));
        acb_pow(w, z, zy, 1000);

        for (rnd = 0; rnd <= ARF_RND_ACCURATE; rnd++)
        {
            acf_pow(r, x, y, 64, rnd);
            check_part("pow", "real", acf_realref(r), acb_realref(w), acb_imagref(w), x, y, 64, rnd);
            check_part("pow", "imag", acf_imagref(r), acb_imagref(w), acb_realref(w), x, y, 64, rnd);
        }

        acf_clear(x);
        acf_clear(y);
        acf_clear(r);
        acb_clear(z);
        acb_clear(zy);
        acb_clear(w);
    }

    /* near the exact points of each function (0, +-1, +-i, the
       multiples of 1/4 and points far out): x = z0 + h with tiny parts of
       h beyond the precision (and in patterns that cancel to first or
       second order), against acb at a precision resolving the distance
       to the rounding boundaries, in the strict modes (the parts lie
       within about |h| or |h|^2 of a dyadic number) */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        acf_t x, r;
        acb_t z, w;
        arf_t t;
        slong prec, k1, k2, wp, which, kind;
        arf_rnd_t rnd;

        acf_init(x);
        acf_init(r);
        acb_init(z);
        acb_init(w);
        arf_init(t);

        prec = 2 + n_randint(state, 100);
        which = n_randint(state, NF);
        k1 = prec + 8 + n_randint(state, 1000);
        k2 = n_randint(state, 3) ? prec + 8 + n_randint(state, 1000) : 2 * k1 + n_randint(state, 3);

        /* h */
        arf_set_si_2exp_si(acf_realref(x), (slong) n_randint(state, 7) - 3, -k1);
        arf_set_si_2exp_si(acf_imagref(x), (slong) n_randint(state, 7) - 3, -k2);
        if (n_randint(state, 2))
            arf_swap(acf_realref(x), acf_imagref(x));

        /* z0 */
        kind = n_randint(state, 6);
        if (kind == 1)
            arf_add_si(acf_realref(x), acf_realref(x), n_randint(state, 2) ? 1 : -1, ARF_PREC_EXACT, ARF_RND_DOWN);
        else if (kind == 2)
            arf_add_si(acf_imagref(x), acf_imagref(x), n_randint(state, 2) ? 1 : -1, ARF_PREC_EXACT, ARF_RND_DOWN);
        else if (kind == 3)
        {
            arf_set_si_2exp_si(t, (slong) n_randint(state, 17) - 8, -2);
            arf_add(acf_realref(x), acf_realref(x), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        }
        else if (kind == 4 && (which == 2 || funcs[which].name[0] == 'a'))
        {
            /* far out (the logarithm and the inverse functions): 2^e (1 +
               h) times a unit */
            slong e = 8 + n_randint(state, 300);
            arf_add_ui(acf_realref(x), acf_realref(x), 1, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(acf_realref(x), acf_realref(x), e);
            arf_mul_2exp_si(acf_imagref(x), acf_imagref(x), e);
            if (n_randint(state, 2))
                arf_swap(acf_realref(x), acf_imagref(x));
            if (n_randint(state, 2))
                arf_neg(acf_realref(x), acf_realref(x));
        }

        if (!arf_is_zero(acf_realref(x)) && !arf_is_zero(acf_imagref(x)))
        {
            wp = 3 * FLINT_MAX(k1, k2) + 4 * prec + 1000;
            arb_set_arf(acb_realref(z), acf_realref(x));
            arb_set_arf(acb_imagref(z), acf_imagref(x));
            funcs[which].g(w, z, wp);
            if (acb_is_finite(w))
            {
                for (rnd = 0; rnd <= ARF_RND_NEAR; rnd++)
                {
                    funcs[which].f(r, x, prec, rnd);
                    check_part(funcs[which].name, "real (near)", acf_realref(r), acb_realref(w), acb_imagref(w), x, NULL, prec, rnd);
                    check_part(funcs[which].name, "imag (near)", acf_imagref(r), acb_imagref(w), acb_realref(w), x, NULL, prec, rnd);
                }
            }
        }

        acf_clear(x);
        acf_clear(r);
        acb_clear(z);
        acb_clear(w);
        arf_clear(t);
    }

    /* pow near its exact cases: x = 2^K i^u (1 + h) with tiny h and y an
       integer, a small dyadic number or complex, and tiny y */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        acf_t x, y, r;
        acb_t z, zy, w;
        slong prec, k1, k2, wp, i;
        arf_rnd_t rnd;

        acf_init(x);
        acf_init(y);
        acf_init(r);
        acb_init(z);
        acb_init(zy);
        acb_init(w);

        prec = 2 + n_randint(state, 100);
        k1 = prec + 40 + n_randint(state, 800);
        k2 = n_randint(state, 3) ? prec + 40 + n_randint(state, 800) : 2 * k1 + n_randint(state, 3);

        if (n_randint(state, 4) != 0)
        {
            arf_set_si_2exp_si(acf_realref(x), (slong) n_randint(state, 7) - 3, -k1);
            arf_set_si_2exp_si(acf_imagref(x), (slong) n_randint(state, 7) - 3, -k2);
            arf_add_ui(acf_realref(x), acf_realref(x), 1, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(acf_realref(x), acf_realref(x), (slong) n_randint(state, 9) - 4);
            arf_mul_2exp_si(acf_imagref(x), acf_imagref(x), (slong) n_randint(state, 9) - 4);
            for (i = n_randint(state, 4); i > 0; i--)
            {
                arf_swap(acf_realref(x), acf_imagref(x));
                arf_neg(acf_realref(x), acf_realref(x));
            }
            switch (n_randint(state, 4))
            {
                case 0: arf_set_si(acf_realref(y), (slong) n_randint(state, 9) - 4); arf_zero(acf_imagref(y)); break;
                case 1: arf_set_si_2exp_si(acf_realref(y), (slong) n_randint(state, 9) - 4, -(slong) n_randint(state, 3)); arf_zero(acf_imagref(y)); break;
                default:
                    arf_set_si_2exp_si(acf_realref(y), (slong) n_randint(state, 9) - 4, -(slong) n_randint(state, 3));
                    arf_set_si_2exp_si(acf_imagref(y), (slong) n_randint(state, 9) - 4, -(slong) n_randint(state, 3));
                    break;
            }
        }
        else
        {
            arf_set_si_2exp_si(acf_realref(x), (slong) n_randint(state, 31) - 15, -(slong) n_randint(state, 4));
            arf_set_si_2exp_si(acf_imagref(x), (slong) n_randint(state, 31) - 15, -(slong) n_randint(state, 4));
            arf_set_si_2exp_si(acf_realref(y), (slong) n_randint(state, 7) - 3, -k1);
            arf_set_si_2exp_si(acf_imagref(y), (slong) n_randint(state, 7) - 3, -k2);
        }

        wp = 3 * FLINT_MAX(k1, k2) + 4 * prec + 1000;
        arb_set_arf(acb_realref(z), acf_realref(x));
        arb_set_arf(acb_imagref(z), acf_imagref(x));
        arb_set_arf(acb_realref(zy), acf_realref(y));
        arb_set_arf(acb_imagref(zy), acf_imagref(y));
        acb_pow(w, z, zy, wp);
        if (acb_is_finite(w) && !(arf_is_zero(acf_realref(x)) && arf_is_zero(acf_imagref(x))))
        {
            for (rnd = 0; rnd <= ARF_RND_NEAR; rnd++)
            {
                acf_pow(r, x, y, prec, rnd);
                check_part("pow", "real (near)", acf_realref(r), acb_realref(w), acb_imagref(w), x, y, prec, rnd);
                check_part("pow", "imag (near)", acf_imagref(r), acb_imagref(w), acb_realref(w), x, y, prec, rnd);
            }
        }

        acf_clear(x);
        acf_clear(y);
        acf_clear(r);
        acb_clear(z);
        acb_clear(zy);
        acb_clear(w);
    }

    /* far beyond any working precision, eps = 2^-(2^62) (a bignum
       exponent): each part just beside a power of two P at 64 bits,
       above (side 1) or below (side -1) it, with the rounding to P or to
       its neighbour P (1 + 2^-63) resp. P (1 - 2^-64):

         exp(eps (1 + i)) = (1 + eps + ..., eps (1 + eps + ...))
         cos(eps (1 + i)) = (1 - eps^4/6 + ..., -eps^2 (1 - eps^4/36))
         log(1 + eps i)  = (eps^2/2 (1 - eps^2/2), eps (1 - eps^2/3))
         (1 + eps i)^2   = (1 - eps^2, 2 eps exactly) */
    {
        acf_t x, y, r;
        arf_t eps, P;
        fmpz_t e;
        int f, part, rnd, st;

        acf_init(x);
        acf_init(y);
        acf_init(r);
        arf_init(eps);
        arf_init(P);
        fmpz_init(e);
        fmpz_one(e);
        fmpz_mul_2exp(e, e, 62);
        fmpz_neg(e, e);
        arf_one(eps);
        arf_mul_2exp_fmpz(eps, eps, e);

        for (f = 0; f < 4; f++)
        {
            /* the powers of two P = 2^c eps^d (signs) and the sides */
            static const int sgn[4][2] = { {1, 1}, {1, -1}, {1, 1}, {1, 1} };
            static const int epsp[4][2] = { {0, 1}, {0, 2}, {2, 1}, {0, 1} };
            static const int c2[4][2] = { {0, 0}, {0, 0}, {-1, 0}, {0, 1} };
            static const int side[4][2] = { {1, 1}, {-1, 1}, {-1, -1}, {-1, 0} };

            if (f == 0 || f == 1)
            {
                arf_set(acf_realref(x), eps);
                arf_set(acf_imagref(x), eps);
            }
            else
            {
                arf_one(acf_realref(x));
                arf_set(acf_imagref(x), eps);
            }

            for (rnd = 0; rnd <= ARF_RND_NEAR; rnd++)
            {
                if (f == 0) st = acf_exp(r, x, 64, rnd);
                else if (f == 1) st = acf_cos(r, x, 64, rnd);
                else if (f == 2) st = acf_log(r, x, 64, rnd);
                else
                {
                    arf_set_si(acf_realref(y), 2);
                    arf_zero(acf_imagref(y));
                    st = acf_pow(r, x, y, 64, rnd);
                }

                for (part = 0; part < 2; part++)
                {
                    int sd = side[f][part], toward;
                    arf_ptr got = part ? acf_imagref(r) : acf_realref(r);

                    arf_set_si(P, sgn[f][part]);
                    if (epsp[f][part] >= 1) arf_mul(P, P, eps, ARF_PREC_EXACT, ARF_RND_DOWN);
                    if (epsp[f][part] >= 2) arf_mul(P, P, eps, ARF_PREC_EXACT, ARF_RND_DOWN);
                    arf_mul_2exp_si(P, P, c2[f][part]);

                    /* the value lies on side sd of P (in value): rounded to
                       the neighbour when the mode moves that way */
                    toward = (sd > 0) ? (rnd == ARF_RND_CEIL || (rnd == ARF_RND_UP && arf_sgn(P) > 0) || (rnd == ARF_RND_DOWN && arf_sgn(P) < 0))
                                      : (rnd == ARF_RND_FLOOR || (rnd == ARF_RND_DOWN && arf_sgn(P) > 0) || (rnd == ARF_RND_UP && arf_sgn(P) < 0));
                    if (sd != 0 && toward)
                    {
                        /* |value| above |P|: |P| (1 + 2^-63); below: |P| (1 - 2^-64) */
                        arf_t T;
                        arf_init(T);
                        if ((sd > 0) == (arf_sgn(P) > 0))
                        {
                            arf_mul_2exp_si(T, P, -63);
                            arf_add(P, P, T, ARF_PREC_EXACT, ARF_RND_DOWN);
                        }
                        else
                        {
                            arf_mul_2exp_si(T, P, -64);
                            arf_sub(P, P, T, ARF_PREC_EXACT, ARF_RND_DOWN);
                        }
                        arf_clear(T);
                    }

                    if (!arf_equal(got, P) || ((st >> part) & 1) != (sd != 0))
                    {
                        flint_printf("FAIL: eps case %d, part %d, rnd = %d, status %d\n", f, part, rnd, st);
                        flint_abort();
                    }
                }
            }
        }

        acf_clear(x);
        acf_clear(y);
        acf_clear(r);
        arf_clear(eps);
        arf_clear(P);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
