/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "decimal.h"
#include "fmpq.h"
#include "gr.h"

/* random finite decfloat with a moderate exponent so that exact rational
   conversion is feasible */
static void
_randtest_moderate(decfloat_t x, flint_rand_t state, gr_ctx_t ctx)
{
    do
    {
        GR_MUST_SUCCEED(decfloat_randtest_special(x, state, ctx));
    }
    while (!DECFLOAT_IS_FINITE(x));

    if (!DECFLOAT_IS_SPECIAL(x))
    {
        slong maxexp = 10 + 100 / DECIMAL_CTX_E(ctx);
        if (n_randint(state, 4) == 0)
            maxexp *= 4;
        fmpz_set_si(&x->exp, (slong) n_randint(state, 2 * maxexp + 1) - maxexp);
    }
}

/* Interval oracle for sqrt: returns 1 if z is the correctly rounded
   square root of q, decided using rational endpoints. */
static int
_check_sqrt(const decfloat_t z, const fmpq_t q, slong prec, int rnd, gr_ctx_t ctx)
{
    /* S = floor(sqrt(q * 10^(2k))) with sticky; sqrt(q) in [S, S+1] / 10^k */
    fmpz_t S, R, t, p;
    fmpq_t lo, hi;
    decfloat_t zl, zh;
    slong k, k0;
    int ok = 0, sticky;

    fmpz_init(S);
    fmpz_init(R);
    fmpz_init(t);
    fmpz_init(p);
    fmpq_init(lo);
    fmpq_init(hi);
    decfloat_init(zl, ctx);
    decfloat_init(zh, ctx);

    /* scale so that S has about prec + 5 digits regardless of the magnitude of q */
    {
        slong E = (slong) fmpz_sizeinbase(fmpq_numref(q), 10) - (slong) fmpz_sizeinbase(fmpq_denref(q), 10);
        k = prec + 5 - E / 2;
        if (k < 1) k = 1;
    }

    for (k0 = k; k < k0 + 1400; k += 40)
    {
        /* t = floor(a * 10^(2k) / b) ; S = isqrt(t) ; sticky if t*b != a*10^(2k) or S^2 != t */
        fmpz_ui_pow_ui(p, 10, 2 * k);
        fmpz_mul(t, fmpq_numref(q), p);
        fmpz_fdiv_qr(t, R, t, fmpq_denref(q));
        sticky = !fmpz_is_zero(R);
        fmpz_sqrtrem(S, R, t);
        sticky |= !fmpz_is_zero(R);

        fmpz_ui_pow_ui(p, 10, k);
        fmpq_set_fmpz_frac(lo, S, p);

        if (!sticky)
        {
            /* exact */
            GR_MUST_SUCCEED(decfloat_set_round_fmpq_reference(zl, lo, prec, rnd, ctx));
            ok = (decfloat_equal(z, zl, ctx) == T_TRUE);
            break;
        }

        fmpz_add_ui(S, S, 1);
        fmpq_set_fmpz_frac(hi, S, p);

        /* sqrt(q) strictly between lo and hi */
        GR_MUST_SUCCEED(decfloat_set_round_fmpq_reference(zl, lo, prec, rnd, ctx));
        GR_MUST_SUCCEED(decfloat_set_round_fmpq_reference(zh, hi, prec, rnd, ctx));

        if (decfloat_equal(zl, zh, ctx) == T_TRUE)
        {
            ok = (decfloat_equal(z, zl, ctx) == T_TRUE);
            break;
        }
    }

    fmpz_clear(S);
    fmpz_clear(R);
    fmpz_clear(t);
    fmpz_clear(p);
    fmpq_clear(lo);
    fmpq_clear(hi);
    decfloat_clear(zl, ctx);
    decfloat_clear(zh, ctx);
    return ok;
}

TEST_FUNCTION_START(decfloat_arith, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t x, y, z, w;
        fmpq_t qx, qy, qz;
        slong prec;
        int rnd, s1, s2, op, aliasing, ok;

        gr_ctx_init_decfloat_randtest(ctx, state, 40);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);

        decfloat_init(x, ctx);
        decfloat_init(y, ctx);
        decfloat_init(z, ctx);
        decfloat_init(w, ctx);
        fmpq_init(qx);
        fmpq_init(qy);
        fmpq_init(qz);

        prec = 1 + n_randint(state, 50);
        if (n_randint(state, 10) == 0)
            prec = DECIMAL_PREC_EXACT;
        rnd = n_randint(state, DECIMAL_RND_NUM);

        _randtest_moderate(x, state, ctx);
        _randtest_moderate(y, state, ctx);

        /* sometimes make y close to x to exercise cancellation */
        if (n_randint(state, 4) == 0)
        {
            decfloat_t t;
            decfloat_init(t, ctx);
            GR_MUST_SUCCEED(decfloat_set_round(y, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
            _randtest_moderate(t, state, ctx);
            if (!DECFLOAT_IS_SPECIAL(y) && !DECFLOAT_IS_SPECIAL(t))
            {
                fmpz_sub_ui(&t->exp, &y->exp, 1 + n_randint(state, 3));
                GR_MUST_SUCCEED(_decfloat_add(y, y, t, n_randint(state, 2), DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, NULL, NULL, ctx));
            }
            decfloat_clear(t, ctx);
        }

        GR_MUST_SUCCEED(decfloat_get_fmpq(qx, x, ctx));
        GR_MUST_SUCCEED(decfloat_get_fmpq(qy, y, ctx));

        op = n_randint(state, 7);
        aliasing = n_randint(state, 4);

        /* z = op(x, y) with random aliasing */
        if (aliasing == 1)
            GR_MUST_SUCCEED(decfloat_set_round(z, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
        else if (aliasing == 2)
            GR_MUST_SUCCEED(decfloat_set_round(z, y, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));

        {
            const decfloat_struct * xa = (aliasing == 1) ? z : x;
            const decfloat_struct * ya = (aliasing == 2) ? z : y;

            if (aliasing == 3)
            {
                ya = xa;
                fmpq_set(qy, qx);
            }

            switch (op)
            {
                case 0:
                    s1 = decfloat_add_round(z, xa, ya, prec, rnd, ctx);
                    fmpq_add(qz, qx, qy);
                    break;
                case 1:
                    s1 = decfloat_sub_round(z, xa, ya, prec, rnd, ctx);
                    fmpq_sub(qz, qx, qy);
                    break;
                case 2:
                    s1 = decfloat_mul_round(z, xa, ya, prec, rnd, ctx);
                    fmpq_mul(qz, qx, qy);
                    break;
                case 3:
                    if (fmpq_is_zero(qy))
                    {
                        s1 = decfloat_zero(z, ctx);
                        fmpq_zero(qz);
                    }
                    else
                    {
                        s1 = decfloat_div_round(z, xa, ya, prec, rnd, ctx);
                        fmpq_div(qz, qx, qy);
                    }
                    break;
                case 4:
                    s1 = decfloat_sqr_round(z, xa, prec, rnd, ctx);
                    fmpq_mul(qz, qx, qx);
                    break;
                case 5:
                    s1 = decfloat_set_round(z, xa, prec, rnd, ctx);
                    fmpq_set(qz, qx);
                    break;
                default:
                    /* sqrt, checked separately below */
                    if (fmpq_sgn(qx) < 0)
                    {
                        fmpq_neg(qx, qx);
                        if (!DECFLOAT_IS_SPECIAL(x)) x->m.size = -x->m.size;
                        if (aliasing == 1 && !DECFLOAT_IS_SPECIAL(z)) z->m.size = -z->m.size;
                    }
                    s1 = decfloat_sqrt_round(z, xa, prec, rnd, ctx);
                    fmpq_set(qz, qx);
                    break;
            }
        }

        if (op == 6)
        {
            if (s1 == GR_SUCCESS)
            {
                if (prec == DECIMAL_PREC_EXACT)
                {
                    /* must be an exact square root */
                    fmpq_t t;
                    fmpq_init(t);
                    GR_MUST_SUCCEED(decfloat_get_fmpq(t, z, ctx));
                    fmpq_mul(t, t, t);
                    ok = fmpq_equal(t, qz);
                    fmpq_clear(t);
                }
                else
                    ok = _check_sqrt(z, qz, prec, rnd, ctx);
            }
            else
            {
                /* failure is only allowed in the exact ring for non-squares */
                ok = (prec == DECIMAL_PREC_EXACT);
                if (ok)
                {
                    fmpq_t t;
                    fmpq_init(t);
                    /* check that q is not a perfect square of a decimal */
                    if (fmpz_is_square(fmpq_numref(qz)) && fmpz_is_square(fmpq_denref(qz)))
                    {
                        fmpz_sqrt(fmpq_numref(t), fmpq_numref(qz));
                        fmpz_sqrt(fmpq_denref(t), fmpq_denref(qz));
                        if (decfloat_set_round_fmpq_reference(w, t, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx) == GR_SUCCESS)
                            ok = 0;
                    }
                    fmpq_clear(t);
                }
            }

            if (!ok)
            {
                flint_printf("FAIL: sqrt\n");
                gr_ctx_println(ctx);
                flint_printf("prec = %wd, rnd = %d, status = %d\n", prec, rnd, s1);
                flint_printf("x = %{gr}\n", x, ctx);
                flint_printf("qx = %{fmpq}\n", qz);
                flint_printf("z = %{gr}\n", z, ctx);
                flint_abort();
            }
        }
        else
        {
            s2 = decfloat_set_round_fmpq_reference(w, qz, prec, rnd, ctx);

            if (s1 != s2 || (s1 == GR_SUCCESS && decfloat_equal(z, w, ctx) != T_TRUE))
            {
                flint_printf("FAIL: op = %d, aliasing = %d\n", op, aliasing);
                gr_ctx_println(ctx);
                flint_printf("prec = %wd, rnd = %d\n", prec, rnd);
                flint_printf("x = %{gr}\n", x, ctx);
                flint_printf("y = %{gr}\n", y, ctx);
                flint_printf("qx = %{fmpq}\n", qx);
                flint_printf("qy = %{fmpq}\n", qy);
                flint_printf("qz = %{fmpq}\n", qz);
                flint_printf("z = %{gr} (%d)\n", z, ctx, s1);
                flint_printf("w = %{gr} (%d)\n", w, ctx, s2);
                flint_abort();
            }
        }

        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        decfloat_clear(z, ctx);
        decfloat_clear(w, ctx);
        fmpq_clear(qx);
        fmpq_clear(qy);
        fmpq_clear(qz);
        gr_ctx_clear(ctx);
    }

    /* large precisions (exercising the FFT / Newton kernels) */
    for (iter = 0; iter < 8 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t x, y, z, w;
        fmpq_t qx, qy, qz;
        slong prec, bits;
        int rnd, s1, s2, op;

        gr_ctx_init_decfloat(ctx, 10, 0);
        prec = 500 + n_randint(state, 6000);
        rnd = n_randint(state, DECIMAL_RND_NUM);

        decfloat_init(x, ctx);
        decfloat_init(y, ctx);
        decfloat_init(z, ctx);
        decfloat_init(w, ctx);
        fmpq_init(qx);
        fmpq_init(qy);
        fmpq_init(qz);

        bits = (slong) (prec * 3.33) + n_randint(state, 100);
        fmpz_randtest(fmpq_numref(qx), state, bits);
        fmpz_ui_pow_ui(fmpq_denref(qx), 10, n_randint(state, prec));
        fmpq_canonicalise(qx);
        fmpz_randtest_not_zero(fmpq_numref(qy), state, bits);
        fmpz_ui_pow_ui(fmpq_denref(qy), 10, n_randint(state, prec));
        fmpq_canonicalise(qy);

        GR_MUST_SUCCEED(decfloat_set_round_fmpq(x, qx, DECIMAL_PREC_EXACT, rnd, ctx));
        GR_MUST_SUCCEED(decfloat_set_round_fmpq(y, qy, DECIMAL_PREC_EXACT, rnd, ctx));

        op = n_randint(state, 4);

        switch (op)
        {
            case 0: s1 = decfloat_mul_round(z, x, y, prec, rnd, ctx); fmpq_mul(qz, qx, qy); break;
            case 1: s1 = decfloat_div_round(z, x, y, prec, rnd, ctx); fmpq_div(qz, qx, qy); break;
            case 2: s1 = decfloat_sqr_round(z, x, prec, rnd, ctx); fmpq_mul(qz, qx, qx); break;
            default:
                fmpq_abs(qx, qx);
                if (!DECFLOAT_IS_SPECIAL(x)) x->m.size = FLINT_ABS(x->m.size);
                s1 = decfloat_sqrt_round(z, x, prec, rnd, ctx);
                fmpq_set(qz, qx);
                break;
        }

        if (op == 3)
        {
            if (s1 != GR_SUCCESS || !_check_sqrt(z, qz, prec, rnd, ctx))
            {
                flint_printf("FAIL: large sqrt, prec = %wd\n", prec);
                flint_abort();
            }
        }
        else
        {
            s2 = decfloat_set_round_fmpq_reference(w, qz, prec, rnd, ctx);
            if (s1 != s2 || (s1 == GR_SUCCESS && decfloat_equal(z, w, ctx) != T_TRUE))
            {
                flint_printf("FAIL: large op = %d, prec = %wd, rnd = %d\n", op, prec, rnd);
                flint_abort();
            }
        }

        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        decfloat_clear(z, ctx);
        decfloat_clear(w, ctx);
        fmpq_clear(qx);
        fmpq_clear(qy);
        fmpq_clear(qz);
        gr_ctx_clear(ctx);
    }

    /* huge exponent differences: results must agree with the eps/sticky logic */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t x, y, z, w;
        slong prec;
        int rnd, s1, s2;

        gr_ctx_init_decfloat_randtest(ctx, state, 40);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        if (DECIMAL_CTX_IS_EXACT(ctx))
            decimal_ctx_set_prec(ctx, 1 + n_randint(state, 40));

        decfloat_init(x, ctx);
        decfloat_init(y, ctx);
        decfloat_init(z, ctx);
        decfloat_init(w, ctx);

        prec = 1 + n_randint(state, 50);
        rnd = n_randint(state, DECIMAL_RND_NUM);

        _randtest_moderate(x, state, ctx);
        _randtest_moderate(y, state, ctx);

        if (!DECFLOAT_IS_SPECIAL(x) && !DECFLOAT_IS_SPECIAL(y))
        {
            /* y very small compared to x */
            fmpz_t t;
            fmpz_init(t);
            fmpz_randtest_unsigned(t, state, 100);
            fmpz_add_ui(t, t, 1000);
            fmpz_sub(&y->exp, &x->exp, t);
            fmpz_clear(t);

            /* x + y should equal x + tiny: compare with x + sign(y) * 10^(-huge) via
               a moderate tiny value computed exactly */
            {
                decfloat_t tiny;
                decfloat_init(tiny, ctx);
                GR_MUST_SUCCEED(decfloat_set_round(tiny, y, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                fmpz_sub_ui(&tiny->exp, &x->exp, FLINT_ABS(x->m.size) + FLINT_ABS(y->m.size) + (prec + DECIMAL_CTX_E(ctx) - 1) / DECIMAL_CTX_E(ctx) + 8);

                s1 = decfloat_add_round(z, x, y, prec, rnd, ctx);
                s2 = decfloat_add_round(w, x, tiny, prec, rnd, ctx);

                if (s1 != s2 || (s1 == GR_SUCCESS && decfloat_equal(z, w, ctx) != T_TRUE))
                {
                    flint_printf("FAIL: add with huge exponent gap\n");
                    gr_ctx_println(ctx);
                    flint_printf("prec = %wd, rnd = %d\n", prec, rnd);
                    flint_printf("x = %{gr}\n", x, ctx);
                    flint_printf("y = %{gr}\n", y, ctx);
                    flint_printf("z = %{gr} (%d)\n", z, ctx, s1);
                    flint_printf("w = %{gr} (%d)\n", w, ctx, s2);
                    flint_abort();
                }

                s1 = decfloat_sub_round(z, x, y, prec, rnd, ctx);
                s2 = decfloat_sub_round(w, x, tiny, prec, rnd, ctx);

                if (s1 != s2 || (s1 == GR_SUCCESS && decfloat_equal(z, w, ctx) != T_TRUE))
                {
                    flint_printf("FAIL: sub with huge exponent gap\n");
                    gr_ctx_println(ctx);
                    flint_printf("prec = %wd, rnd = %d\n", prec, rnd);
                    flint_printf("x = %{gr}\n", x, ctx);
                    flint_printf("y = %{gr}\n", y, ctx);
                    flint_printf("z = %{gr} (%d)\n", z, ctx, s1);
                    flint_printf("w = %{gr} (%d)\n", w, ctx, s2);
                    flint_abort();
                }

                decfloat_clear(tiny, ctx);
            }
        }

        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        decfloat_clear(z, ctx);
        decfloat_clear(w, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
