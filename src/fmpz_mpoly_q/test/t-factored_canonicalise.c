/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "mpoly.h"
#include "fmpz_mpoly_q.h"

#ifndef FMPZ_MPOLY_Q_FACTORED_TEST_HELPERS
#define FMPZ_MPOLY_Q_FACTORED_TEST_HELPERS

static void
_factored_randtest(fmpz_mpoly_q_factored_t A, flint_rand_t state, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_factored_randtest(A, state, 1 + n_randint(state, 4),
        1 + n_randint(state, 20), 1 + n_randint(state, 3), ctx);
}

/* Generate B, which is likely to share factors with A. */
static void
_factored_randtest_related(fmpz_mpoly_q_factored_t B, const fmpz_mpoly_q_factored_t A,
    flint_rand_t state, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_t a, b, c;
    slong k;

    _factored_randtest(B, state, ctx);

    if (n_randint(state, 4) == 0 || fmpz_mpoly_q_factored_is_zero(A, ctx))
        return;

    fmpz_mpoly_q_init(a, ctx);
    fmpz_mpoly_q_init(b, ctx);
    fmpz_mpoly_q_init(c, ctx);

    fmpz_mpoly_q_factored_get_fmpz_mpoly_q(a, A, ctx);
    fmpz_mpoly_q_factored_get_fmpz_mpoly_q(c, B, ctx);

    switch (n_randint(state, 4))
    {
        case 0:
            /* b = c + a */
            fmpz_mpoly_q_add(b, c, a, ctx);
            break;
        case 1:
            /* b = c / den(a)^k */
            fmpz_mpoly_set(fmpz_mpoly_q_numref(b), fmpz_mpoly_q_denref(a), ctx);
            fmpz_mpoly_one(fmpz_mpoly_q_denref(b), ctx);
            for (k = n_randint(state, 2); k >= 0; k--)
                fmpz_mpoly_q_div(c, c, b, ctx);
            fmpz_mpoly_q_swap(b, c, ctx);
            break;
        case 2:
            /* b = c * a^(+/-1) */
            if (n_randint(state, 2))
                fmpz_mpoly_q_div(b, c, a, ctx);
            else
                fmpz_mpoly_q_mul(b, c, a, ctx);
            break;
        default:
            /* b = a * (a + c) */
            fmpz_mpoly_q_add(b, a, c, ctx);
            fmpz_mpoly_q_mul(b, b, a, ctx);
            break;
    }

    fmpz_mpoly_q_factored_set_fmpz_mpoly_q(B, b, ctx);

    fmpz_mpoly_q_clear(a, ctx);
    fmpz_mpoly_q_clear(b, ctx);
    fmpz_mpoly_q_clear(c, ctx);
}

static void
_factored_check(const fmpz_mpoly_q_factored_t C, const fmpz_mpoly_q_t ref,
    const fmpz_mpoly_q_factored_t A, const fmpz_mpoly_q_factored_t B,
    const char * op, int check_canonical, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_t c;
    fmpz_mpoly_q_factored_t D;
    int ok;

    fmpz_mpoly_q_init(c, ctx);
    fmpz_mpoly_q_factored_init(D, ctx);

    fmpz_mpoly_q_factored_get_fmpz_mpoly_q(c, C, ctx);
    fmpz_mpoly_q_factored_set_fmpz_mpoly_q(D, ref, ctx);

    ok = fmpz_mpoly_q_equal(c, ref, ctx) && fmpz_mpoly_q_factored_equal(C, D, ctx);

    if (ok && check_canonical)
        ok = fmpz_mpoly_q_factored_is_canonical(C, ctx);

    if (!ok)
    {
        flint_printf("FAIL: %s\n", op);
        flint_printf("A = "); fmpz_mpoly_q_factored_print_pretty(A, NULL, ctx); flint_printf("\n\n");
        if (B != NULL)
        {
            flint_printf("B = "); fmpz_mpoly_q_factored_print_pretty(B, NULL, ctx); flint_printf("\n\n");
        }
        flint_printf("C = "); fmpz_mpoly_q_factored_print_pretty(C, NULL, ctx); flint_printf("\n\n");
        flint_printf("D = "); fmpz_mpoly_q_factored_print_pretty(D, NULL, ctx); flint_printf("\n\n");
        flint_printf("ref = "); fmpz_mpoly_q_print_pretty(ref, NULL, ctx); flint_printf("\n\n");
        flint_abort();
    }

    fmpz_mpoly_q_clear(c, ctx);
    fmpz_mpoly_q_factored_clear(D, ctx);
}

#endif

TEST_FUNCTION_START(fmpz_mpoly_q_factored_canonicalise, state)
{
    slong iter;

    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_factored_t A, B;
        fmpz_mpoly_q_t a, b;
        fmpz_mpoly_t t;
        fmpz_t c;
        char * s;

        fmpz_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), mpoly_ordering_randtest(state));

        fmpz_mpoly_q_factored_init(A, ctx);
        fmpz_mpoly_q_factored_init(B, ctx);
        fmpz_mpoly_q_init(a, ctx);
        fmpz_mpoly_q_init(b, ctx);
        fmpz_mpoly_init(t, ctx);
        fmpz_init(c);

        _factored_randtest(A, state, ctx);

        if (!fmpz_mpoly_q_factored_is_canonical(A, ctx))
        {
            flint_printf("FAIL: randtest not canonical\n");
            flint_printf("A = "); fmpz_mpoly_q_factored_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_abort();
        }

        /* round trip via fmpz_mpoly_q */
        fmpz_mpoly_q_factored_get_fmpz_mpoly_q(a, A, ctx);

        if (!fmpz_mpoly_q_is_canonical(a, ctx))
        {
            flint_printf("FAIL: get_fmpz_mpoly_q not canonical\n");
            flint_printf("A = "); fmpz_mpoly_q_factored_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_abort();
        }

        _factored_randtest(B, state, ctx);
        fmpz_mpoly_q_factored_set_fmpz_mpoly_q(B, a, ctx);

        if (!fmpz_mpoly_q_factored_equal(A, B, ctx))
        {
            flint_printf("FAIL: round trip\n");
            flint_printf("A = "); fmpz_mpoly_q_factored_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mpoly_q_factored_print_pretty(B, NULL, ctx); flint_printf("\n\n");
            flint_abort();
        }

        /* multiply numerator and denominator by a common random
           polynomial and constant (stored with a random sign, unfactored),
           then canonicalise */
        fmpz_mpoly_q_factored_set(B, A, ctx);
        do {
            fmpz_mpoly_randtest_bound(t, state, 1 + n_randint(state, 4), 1 + n_randint(state, 20), 1 + n_randint(state, 3), ctx);
        } while (fmpz_mpoly_is_zero(t, ctx));
        fmpz_randtest_not_zero(c, state, 20);

        fmpz_mpoly_mul(fmpz_mpoly_q_factored_numref(B), fmpz_mpoly_q_factored_numref(B), t, ctx);
        fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_q_factored_numref(B), fmpz_mpoly_q_factored_numref(B), c, ctx);
        fmpz_mpoly_factor_append_ui(fmpz_mpoly_q_factored_denref(B), t, 1, ctx);
        fmpz_mul(B->den.constant, B->den.constant, c);
        if (n_randint(state, 2))
        {
            fmpz_mpoly_neg(fmpz_mpoly_q_factored_numref(B), fmpz_mpoly_q_factored_numref(B), ctx);
            fmpz_neg(B->den.constant, B->den.constant);
        }
        if (n_randint(state, 2))
        {
            /* a factor with exponent zero should be ignored */
            fmpz_mpoly_factor_append_ui(fmpz_mpoly_q_factored_denref(B), t, 0, ctx);
        }

        fmpz_mpoly_q_factored_canonicalise(B, ctx);

        if (!fmpz_mpoly_q_factored_equal(A, B, ctx) || !fmpz_mpoly_q_factored_is_canonical(B, ctx))
        {
            flint_printf("FAIL: canonicalise\n");
            flint_printf("A = "); fmpz_mpoly_q_factored_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mpoly_q_factored_print_pretty(B, NULL, ctx); flint_printf("\n\n");
            flint_abort();
        }

        /* negative exponents and constant_den: B = A * t^k * c, k >= 1 */
        {
            fmpz_mpoly_q_factored_t E, T;
            slong k = 1 + n_randint(state, 2);

            fmpz_mpoly_q_factored_init(E, ctx);
            fmpz_mpoly_q_factored_init(T, ctx);

            fmpz_mpoly_q_factored_set_fmpz_mpoly(T, t, ctx);
            fmpz_mpoly_q_factored_pow_si(T, T, k, ctx);
            fmpz_mpoly_q_factored_mul(E, A, T, ctx);
            fmpz_mpoly_q_factored_set_fmpz(T, c, ctx);
            fmpz_mpoly_q_factored_mul(E, E, T, ctx);

            fmpz_mpoly_q_factored_set(B, A, ctx);
            fmpz_mpoly_factor_append_ui(fmpz_mpoly_q_factored_denref(B), t, 0, ctx);
            fmpz_set_si(B->den.exp + B->den.num - 1, -k);
            fmpz_set(B->den.constant_den, c);
            fmpz_mpoly_q_factored_canonicalise(B, ctx);

            if (!fmpz_mpoly_q_factored_equal(B, E, ctx) || !fmpz_mpoly_q_factored_is_canonical(B, ctx))
            {
                flint_printf("FAIL: canonicalise (negative exponents)\n");
                flint_printf("A = "); fmpz_mpoly_q_factored_print_pretty(A, NULL, ctx); flint_printf("\n\n");
                flint_printf("B = "); fmpz_mpoly_q_factored_print_pretty(B, NULL, ctx); flint_printf("\n\n");
                flint_printf("E = "); fmpz_mpoly_q_factored_print_pretty(E, NULL, ctx); flint_printf("\n\n");
                flint_abort();
            }

            fmpz_mpoly_q_factored_clear(E, ctx);
            fmpz_mpoly_q_factored_clear(T, ctx);
        }

        /* conversions from simple types agree with fmpz_mpoly_q */
        {
            fmpq_t q;
            slong v;

            fmpq_init(q);
            fmpq_randtest(q, state, 20);
            v = n_randint(state, fmpz_mpoly_ctx_nvars(ctx));

            switch (n_randint(state, 6))
            {
                case 0:
                    fmpz_mpoly_q_set_fmpq(b, q, ctx);
                    fmpz_mpoly_q_factored_set_fmpq(B, q, ctx);
                    break;
                case 1:
                    fmpz_mpoly_q_set_fmpz(b, fmpq_numref(q), ctx);
                    fmpz_mpoly_q_factored_set_fmpz(B, fmpq_numref(q), ctx);
                    break;
                case 2:
                    v = (slong) n_randint(state, 5) - 2;
                    fmpz_mpoly_q_set_si(b, v, ctx);
                    fmpz_mpoly_q_factored_set_si(B, v, ctx);
                    break;
                case 3:
                    fmpz_mpoly_q_gen(b, v, ctx);
                    fmpz_mpoly_q_factored_gen(B, v, ctx);
                    break;
                case 4:
                    fmpz_mpoly_q_one(b, ctx);
                    fmpz_mpoly_q_factored_one(B, ctx);
                    break;
                default:
                    fmpz_mpoly_q_one(b, ctx);
                    fmpz_mpoly_set(fmpz_mpoly_q_numref(b), t, ctx);
                    fmpz_mpoly_q_factored_set_fmpz_mpoly(B, t, ctx);
                    break;
            }

            fmpz_mpoly_q_factored_get_fmpz_mpoly_q(a, B, ctx);

            if (!fmpz_mpoly_q_equal(a, b, ctx) || !fmpz_mpoly_q_factored_is_canonical(B, ctx) ||
                fmpz_mpoly_q_is_zero(b, ctx) != fmpz_mpoly_q_factored_is_zero(B, ctx) ||
                fmpz_mpoly_q_is_one(b, ctx) != fmpz_mpoly_q_factored_is_one(B, ctx))
            {
                flint_printf("FAIL: conversions\n");
                flint_printf("B = "); fmpz_mpoly_q_factored_print_pretty(B, NULL, ctx); flint_printf("\n\n");
                flint_printf("b = "); fmpz_mpoly_q_print_pretty(b, NULL, ctx); flint_printf("\n\n");
                flint_abort();
            }

            fmpz_mpoly_q_factored_get_fmpz_mpoly_q(a, A, ctx);
            fmpq_clear(q);
        }

        /* the string representation should parse to the same value */
        s = fmpz_mpoly_q_factored_get_str_pretty(A, NULL, ctx);
        {
            const char * vars[] = { "x1", "x2", "x3", "x4" };
            if (fmpz_mpoly_q_set_str_pretty(b, s, vars, ctx) != 0 || !fmpz_mpoly_q_equal(a, b, ctx))
            {
                flint_printf("FAIL: get_str_pretty\n");
                flint_printf("s = %s\n\n", s);
                flint_printf("a = "); fmpz_mpoly_q_print_pretty(a, NULL, ctx); flint_printf("\n\n");
                flint_abort();
            }
        }
        flint_free(s);

        fmpz_mpoly_q_factored_clear(A, ctx);
        fmpz_mpoly_q_factored_clear(B, ctx);
        fmpz_mpoly_q_clear(a, ctx);
        fmpz_mpoly_q_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_clear(c);

        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
