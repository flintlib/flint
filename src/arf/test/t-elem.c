/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "elem_check.h"

/* The elementary functions of arf against balls of arb (or of the gr
   context of arb, for the functions arb lacks): correctly rounded in the
   strict modes, within the tolerance of the relaxed modes, with aliasing
   (elem_check.h); the functions in turn.  Then the sine-cosine pairs
   against the separate functions. */

/* 2^x: exact for integers, else exp(x log 2) */
static void
_arb_exp2(arb_ptr res, const arb_t x, slong prec)
{
    if (arf_is_finite(arb_midref(x)) && arf_cmpabs_2exp_si(arb_midref(x), 100000) >= 0)
    {
        arb_zero_pm_inf(res);
    }
    else if (arb_is_exact(x) && arf_is_int(arb_midref(x)) && arf_is_finite(arb_midref(x)))
    {
        fmpz_t k;
        fmpz_init(k);
        arf_get_fmpz(k, arb_midref(x), ARF_RND_DOWN);
        arb_one(res);
        arb_mul_2exp_fmpz(res, res, k);
        fmpz_clear(k);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_const_log2(t, prec + 10);
        arb_mul(t, t, x, prec + 10);
        arb_exp(res, t, prec);
        arb_clear(t);
    }
}

/* log2(x) = log(x) / log 2, exact for powers of two */
static void
_arb_log2(arb_ptr res, const arb_t x, slong prec)
{
    if (arb_is_exact(x) && arf_is_finite(arb_midref(x)) && arf_sgn(arb_midref(x)) > 0
        && (ARF_IS_POW2(arb_midref(x))))
    {
        fmpz_t k;
        fmpz_init(k);
        fmpz_sub_ui(k, ARF_EXPREF(arb_midref(x)), 1);
        arb_set_round_fmpz(res, k, prec);
        fmpz_clear(k);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_const_log2(t, prec + 10);
        arb_log(res, x, prec + 10);
        arb_div(res, res, t, prec);
        arb_clear(t);
    }
}

/* the functions: the reference as an arb function, or (NULL) as a method
   of the gr context of arb; trig: the magnitude limit (0 of exp, 1 of the
   trigonometric functions, 2 of exp on x log 10) */
static const struct
{
    const char * name;
    arf_elem_func_t f;
    arb_elem_func_t g;
    int method;
    int trig;
}
arf_elem_funcs[] =
{
    { "arf_exp", arf_exp, arb_exp, 0, 0 },
    { "arf_expm1", arf_expm1, arb_expm1, 0, 0 },
    { "arf_exp2", arf_exp2, _arb_exp2, 0, 0 },
    { "arf_exp10", arf_exp10, NULL, GR_METHOD_EXP10, 2 },
    { "arf_log", arf_log, arb_log, 0, 0 },
    { "arf_log1p", arf_log1p, arb_log1p, 0, 0 },
    { "arf_log2", arf_log2, _arb_log2, 0, 0 },
    { "arf_log10", arf_log10, NULL, GR_METHOD_LOG10, 0 },
    { "arf_sin", arf_sin, arb_sin, 0, 1 },
    { "arf_cos", arf_cos, arb_cos, 0, 1 },
    { "arf_tan", arf_tan, arb_tan, 0, 1 },
    { "arf_cot", arf_cot, NULL, GR_METHOD_COT, 1 },
    { "arf_sec", arf_sec, NULL, GR_METHOD_SEC, 1 },
    { "arf_csc", arf_csc, NULL, GR_METHOD_CSC, 1 },
    { "arf_sinc", arf_sinc, NULL, GR_METHOD_SINC, 1 },
    { "arf_sin_pi", arf_sin_pi, arb_sin_pi, 0, 1 },
    { "arf_cos_pi", arf_cos_pi, arb_cos_pi, 0, 1 },
    { "arf_tan_pi", arf_tan_pi, arb_tan_pi, 0, 1 },
    { "arf_cot_pi", arf_cot_pi, NULL, GR_METHOD_COT_PI, 1 },
    { "arf_sec_pi", arf_sec_pi, NULL, GR_METHOD_SEC_PI, 1 },
    { "arf_csc_pi", arf_csc_pi, NULL, GR_METHOD_CSC_PI, 1 },
    { "arf_sinc_pi", arf_sinc_pi, NULL, GR_METHOD_SINC_PI, 0 },
    { "arf_sinh", arf_sinh, arb_sinh, 0, 0 },
    { "arf_cosh", arf_cosh, arb_cosh, 0, 0 },
    { "arf_tanh", arf_tanh, arb_tanh, 0, 0 },
    { "arf_coth", arf_coth, arb_coth, 0, 0 },
    { "arf_sech", arf_sech, arb_sech, 0, 0 },
    { "arf_csch", arf_csch, arb_csch, 0, 0 },
    { "arf_atan", arf_atan, arb_atan, 0, 0 },
    { "arf_asin", arf_asin, NULL, GR_METHOD_ASIN, 0 },
    { "arf_acos", arf_acos, NULL, GR_METHOD_ACOS, 0 },
    { "arf_acot", arf_acot, NULL, GR_METHOD_ACOT, 0 },
    { "arf_asec", arf_asec, NULL, GR_METHOD_ASEC, 0 },
    { "arf_acsc", arf_acsc, NULL, GR_METHOD_ACSC, 0 },
    { "arf_atan_pi", arf_atan_pi, NULL, GR_METHOD_ATAN_PI, 0 },
    { "arf_asin_pi", arf_asin_pi, NULL, GR_METHOD_ASIN_PI, 0 },
    { "arf_acos_pi", arf_acos_pi, NULL, GR_METHOD_ACOS_PI, 0 },
    { "arf_asinh", arf_asinh, NULL, GR_METHOD_ASINH, 0 },
    { "arf_acosh", arf_acosh, NULL, GR_METHOD_ACOSH, 0 },
    { "arf_atanh", arf_atanh, NULL, GR_METHOD_ATANH, 0 },
    { "arf_acoth", arf_acoth, NULL, GR_METHOD_ACOTH, 0 },
    { "arf_asech", arf_asech, NULL, GR_METHOD_ASECH, 0 },
    { "arf_acsch", arf_acsch, NULL, GR_METHOD_ACSCH, 0 },
};

#define ARF_ELEM_NFUNCS ((slong) (sizeof(arf_elem_funcs) / sizeof(arf_elem_funcs[0])))

typedef int (* arf_elem_pair_func_t)(arf_ptr, arf_ptr, arf_srcptr, slong, arf_rnd_t);
typedef void (* arb_elem_pair_func_t)(arb_ptr, arb_ptr, const arb_t, slong);

/* the pair f(x) = (s, c) against the separate functions fs, fc (equal
   results in the strict modes; in the relaxed modes, which may round the
   pair and the separate values differently, against the reference g), and
   with aliasing */
static void
arf_elem_check_pair(const char * name, arf_elem_pair_func_t f, arf_elem_func_t fs, arf_elem_func_t fc,
    arb_elem_pair_func_t g, slong maglim_trig, flint_rand_t state)
{
    arf_t x, s, c, s2, c2;
    slong prec, maglim;
    arf_rnd_t rnd;
    int r1, r2;

    arf_init(x);
    arf_init(s);
    arf_init(c);
    arf_init(s2);
    arf_init(c2);

    prec = 2 + n_randint(state, 300);
    if (n_randint(state, 8) == 0)
        prec = 2 + n_randint(state, 3000);
    rnd = (arf_rnd_t) n_randint(state, 7);
    maglim = maglim_trig ? FLINT_MAX(65536, 4 * prec) : FLINT_MAX(128, 2 * prec);
    arf_elem_randtest_arg(x, state, prec);

    r1 = f(s, c, x, prec, rnd);
    r2 = fs(s2, x, prec, rnd);
    fc(c2, x, prec, rnd);

    if (rnd >= ARF_RND_FAST)
    {
        arb_t bs, bc, bx;
        arb_init(bs);
        arb_init(bc);
        arb_init(bx);
        arb_set_arf(bx, x);
        g(bs, bc, bx, prec + 64);
        arf_elem_verify(name, x, s, r1, bs, prec, rnd, maglim);
        arf_elem_verify(name, x, c, r1, bc, prec, rnd, maglim);
        arb_clear(bs);
        arb_clear(bc);
        arb_clear(bx);
        arf_set(s2, s);
        arf_set(c2, c);
    }

    if (!arf_equal(s, s2) || !arf_equal(c, c2) || (r1 != r2 && !arf_is_nan(s) && rnd < ARF_RND_FAST))
    {
        flint_printf("FAIL: %s, prec = %wd, rnd = %d\n", name, prec, rnd);
        flint_printf("x = "); arf_debug(x); flint_printf("\n");
        flint_printf("s = "); arf_debug(s); flint_printf("\n");
        flint_printf("s2 = "); arf_debug(s2); flint_printf("\n");
        flint_printf("c = "); arf_debug(c); flint_printf("\n");
        flint_printf("c2 = "); arf_debug(c2); flint_printf("\n");
        flint_abort();
    }

    /* aliasing */
    arf_set(s, x);
    f(s, c, s, prec, rnd);
    if (!arf_equal(s, s2) || !arf_equal(c, c2))
    {
        flint_printf("FAIL (aliasing): %s, prec = %wd, rnd = %d\n", name, prec, rnd);
        flint_printf("x = "); arf_debug(x); flint_printf("\n");
        flint_abort();
    }

    arf_clear(x);
    arf_clear(s);
    arf_clear(c);
    arf_clear(s2);
    arf_clear(c2);
}

TEST_FUNCTION_START(arf_elem, state)
{
    slong iter, i, iters = 2000 * flint_test_multiplier();

    for (i = 0; i < ARF_ELEM_NFUNCS; i++)
    {
        if (arf_elem_funcs[i].g != NULL)
            arf_elem_test(arf_elem_funcs[i].name, arf_elem_funcs[i].f, arf_elem_funcs[i].g,
                arf_elem_funcs[i].trig, state, iters);
        else
            arf_elem_test_gr(arf_elem_funcs[i].name, arf_elem_funcs[i].f, arf_elem_funcs[i].method,
                arf_elem_funcs[i].trig, state, iters);
    }

    for (iter = 0; iter < iters; iter++)
    {
        arf_elem_check_pair("arf_sin_cos", arf_sin_cos, arf_sin, arf_cos, arb_sin_cos, 1, state);
        arf_elem_check_pair("arf_sin_cos_pi", arf_sin_cos_pi, arf_sin_pi, arf_cos_pi, arb_sin_cos_pi, 1, state);
        arf_elem_check_pair("arf_sinh_cosh", arf_sinh_cosh, arf_sinh, arf_cosh, arb_sinh_cosh, 0, state);
    }

    /* the exact integer values at every precision: log10(10^k) = k,
       log2(2^k) = k, exp10(k) = 10^k, exp2(k) = 2^k, rounded to prec
       (the value of the relaxed modes rounded as their strict mode, as
       these exact values are) */
    {
        arf_t x, r, t;
        fmpz_t z;
        slong k, prec;
        int rnd, s1, s2, f;

        arf_init(x);
        arf_init(r);
        arf_init(t);
        fmpz_init(z);

        for (k = 0; k <= 40; k++)
        for (prec = 1; prec <= 8; prec++)
        for (rnd = 0; rnd <= ARF_RND_ACCURATE; rnd++)
        for (f = 0; f < 4; f++)
        {
            arf_rnd_t sr = arf_rnd_relaxed_to_strict((arf_rnd_t) rnd);

            if (f == 0 || f == 2)
            {
                fmpz_ui_pow_ui(z, 10, k);
                arf_set_fmpz(x, z);
            }
            else
                arf_set_ui_2exp_si(x, 1, k);

            if (f == 0)
                s1 = arf_log10(r, x, prec, (arf_rnd_t) rnd);
            else if (f == 1)
                s1 = arf_log2(r, x, prec, (arf_rnd_t) rnd);

            if (f <= 1)
                s2 = arf_set_round_si(t, k, prec, sr);
            else
            {
                arf_set_si(x, k);
                s1 = (f == 2) ? arf_exp10(r, x, prec, (arf_rnd_t) rnd) : arf_exp2(r, x, prec, (arf_rnd_t) rnd);
                if (f == 2)
                    s2 = arf_set_round_fmpz(t, z, prec, sr);
                else
                    s2 = arf_set_round(t, (arf_set_ui_2exp_si(t, 1, k), t), prec, sr);
            }

            /* (in the relaxed modes, arf_pow and so arf_exp10 do not look
               for exact powers: the flag may be 1 for an exact value) */
            if (!arf_equal(r, t) || (s1 == 0 && s2 != 0)
                || (rnd < ARF_RND_FAST && s1 != 0 && s2 == 0))
            {
                flint_printf("FAIL (exact): f = %d, k = %wd, prec = %wd, rnd = %d\n", f, k, prec, rnd);
                flint_printf("r = "); arf_printd(r, 20); flint_printf(" %d\n", s1);
                flint_printf("t = "); arf_printd(t, 20); flint_printf(" %d\n", s2);
                flint_abort();
            }
        }

        arf_clear(x);
        arf_clear(r);
        arf_clear(t);
        fmpz_clear(z);
    }

    TEST_FUNCTION_END(state);
}
