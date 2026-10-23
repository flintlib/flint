/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
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

TEST_FUNCTION_START(fmpz_mpoly_q_factored_scalar, state)
{
    slong iter;

    for (iter = 0; iter < 1500 * flint_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_factored_t A, C;
        fmpz_mpoly_q_t a, b, c;
        fmpz_mpoly_t p;
        fmpq_t q;
        slong si;
        int op, kind, alias;
        char name[64];

        fmpz_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), mpoly_ordering_randtest(state));

        fmpz_mpoly_q_factored_init(A, ctx);
        fmpz_mpoly_q_factored_init(C, ctx);
        fmpz_mpoly_q_init(a, ctx);
        fmpz_mpoly_q_init(b, ctx);
        fmpz_mpoly_q_init(c, ctx);
        fmpz_mpoly_init(p, ctx);
        fmpq_init(q);

        _factored_randtest(A, state, ctx);
        fmpz_mpoly_q_factored_get_fmpz_mpoly_q(a, A, ctx);

        op = n_randint(state, 4);     /* add, sub, mul, div */
        kind = n_randint(state, 4);   /* fmpz, fmpq, si, fmpz_mpoly */
        alias = n_randint(state, 2);

        /* random scalar; sometimes share factors with A's denominator */
        fmpq_randtest(q, state, 1 + n_randint(state, 40));
        si = z_randtest(state);
        if (n_randint(state, 4) == 0)
            fmpz_mpoly_set(p, fmpz_mpoly_q_denref(a), ctx);
        else
            fmpz_mpoly_randtest_bound(p, state, 1 + n_randint(state, 4), 1 + n_randint(state, 20), 1 + n_randint(state, 3), ctx);
        if (n_randint(state, 2))
            fmpz_mpoly_mul(p, p, fmpz_mpoly_q_denref(a), ctx);

        if (kind == 0)
            fmpz_one(fmpq_denref(q));

        if (op == 3)
        {
            if (fmpq_is_zero(q)) fmpq_one(q);
            if (si == 0) si = 1;
            if (fmpz_mpoly_is_zero(p, ctx)) fmpz_mpoly_one(p, ctx);
        }

        /* reference b = scalar as fmpz_mpoly_q */
        if (kind == 0 || kind == 1)
            fmpz_mpoly_q_set_fmpq(b, q, ctx);
        else if (kind == 2)
            fmpz_mpoly_q_set_si(b, si, ctx);
        else
        {
            fmpz_mpoly_set(fmpz_mpoly_q_numref(b), p, ctx);
            fmpz_mpoly_one(fmpz_mpoly_q_denref(b), ctx);
        }

        if (op == 0) fmpz_mpoly_q_add(c, a, b, ctx);
        else if (op == 1) fmpz_mpoly_q_sub(c, a, b, ctx);
        else if (op == 2) fmpz_mpoly_q_mul(c, a, b, ctx);
        else fmpz_mpoly_q_div(c, a, b, ctx);

        if (alias)
            fmpz_mpoly_q_factored_set(C, A, ctx);
        else
            _factored_randtest(C, state, ctx);

#define FMPZ_MPOLY_Q_FACTORED_DISPATCH(suffix, arg) \
        if (op == 0) fmpz_mpoly_q_factored_add_ ## suffix(C, alias ? C : A, arg, ctx); \
        else if (op == 1) fmpz_mpoly_q_factored_sub_ ## suffix(C, alias ? C : A, arg, ctx); \
        else if (op == 2) fmpz_mpoly_q_factored_mul_ ## suffix(C, alias ? C : A, arg, ctx); \
        else fmpz_mpoly_q_factored_div_ ## suffix(C, alias ? C : A, arg, ctx);

        if (kind == 0)
        {
            FMPZ_MPOLY_Q_FACTORED_DISPATCH(fmpz, fmpq_numref(q))
        }
        else if (kind == 1)
        {
            FMPZ_MPOLY_Q_FACTORED_DISPATCH(fmpq, q)
        }
        else if (kind == 2)
        {
            FMPZ_MPOLY_Q_FACTORED_DISPATCH(si, si)
        }
        else
        {
            FMPZ_MPOLY_Q_FACTORED_DISPATCH(fmpz_mpoly, p)
        }

#undef FMPZ_MPOLY_Q_FACTORED_DISPATCH

        flint_sprintf(name, "scalar op = %d, kind = %d, alias = %d", op, kind, alias);
        _factored_check(C, c, A, NULL, name, iter % 4 == 0, ctx);

        fmpz_mpoly_q_factored_clear(A, ctx);
        fmpz_mpoly_q_factored_clear(C, ctx);
        fmpz_mpoly_q_clear(a, ctx);
        fmpz_mpoly_q_clear(b, ctx);
        fmpz_mpoly_q_clear(c, ctx);
        fmpz_mpoly_clear(p, ctx);
        fmpq_clear(q);

        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
