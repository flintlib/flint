/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
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

TEST_FUNCTION_START(fmpz_mpoly_q_factored_get_set_str, state)
{
    slong iter;

    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_factored_t A, B;
        fmpz_mpoly_q_t a, b;
        const char * vars1[] = { "x1", "x2", "x3", "x4" };
        const char * vars2[] = { "x", "y", "zeta", "w_2" };
        const char ** vars;
        char * s;
        int ret;

        fmpz_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), mpoly_ordering_randtest(state));

        fmpz_mpoly_q_factored_init(A, ctx);
        fmpz_mpoly_q_factored_init(B, ctx);
        fmpz_mpoly_q_init(a, ctx);
        fmpz_mpoly_q_init(b, ctx);

        vars = n_randint(state, 2) ? vars1 : vars2;

        _factored_randtest(A, state, ctx);
        _factored_randtest(B, state, ctx);

        /* round trip in factored form */
        s = fmpz_mpoly_q_factored_get_str_pretty(A, (n_randint(state, 2) || vars == vars2) ? vars : NULL, ctx);
        ret = fmpz_mpoly_q_factored_set_str_pretty(B, s, vars, ctx);

        if (ret != 0 || !fmpz_mpoly_q_factored_equal(A, B, ctx) || !fmpz_mpoly_q_factored_is_canonical(B, ctx))
        {
            flint_printf("FAIL: round trip\n");
            flint_printf("s = %s\n", s);
            flint_printf("A = "); fmpz_mpoly_q_factored_print_pretty(A, vars, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mpoly_q_factored_print_pretty(B, vars, ctx); flint_printf("\n\n");
            flint_abort();
        }

        flint_free(s);

        /* parsing an expanded fmpz_mpoly_q string */
        fmpz_mpoly_q_factored_get_fmpz_mpoly_q(a, A, ctx);
        s = fmpz_mpoly_q_get_str_pretty(a, vars, ctx);
        _factored_randtest(B, state, ctx);
        ret = fmpz_mpoly_q_factored_set_str_pretty(B, s, vars, ctx);

        if (ret != 0 || !fmpz_mpoly_q_factored_equal(A, B, ctx))
        {
            flint_printf("FAIL: parse expanded\n");
            flint_printf("s = %s\n", s);
            flint_printf("A = "); fmpz_mpoly_q_factored_print_pretty(A, vars, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mpoly_q_factored_print_pretty(B, vars, ctx); flint_printf("\n\n");
            flint_abort();
        }

        flint_free(s);

        fmpz_mpoly_q_factored_clear(A, ctx);
        fmpz_mpoly_q_factored_clear(B, ctx);
        fmpz_mpoly_q_clear(a, ctx);
        fmpz_mpoly_q_clear(b, ctx);

        fmpz_mpoly_ctx_clear(ctx);
    }

    /* fixed examples */
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_factored_t A;
        const char * vars[] = { "x", "y" };
        char * s;

        fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
        fmpz_mpoly_q_factored_init(A, ctx);

        if (fmpz_mpoly_q_factored_set_str_pretty(A, "(x^2-1)/(6*(x+1)*(2*y-2)^2) + 1/(3*(y-1))", vars, ctx) != 0)
        {
            flint_printf("FAIL: parse example\n");
            flint_abort();
        }

        s = fmpz_mpoly_q_factored_get_str_pretty(A, vars, ctx);
        if (strcmp(s, "(x+8*y-9)/(24*(y-1)^2)") != 0)
        {
            flint_printf("FAIL: example\n");
            flint_printf("s = %s\n", s);
            flint_abort();
        }
        flint_free(s);

        if (fmpz_mpoly_q_factored_set_str_pretty(A, "x/(y-y)", vars, ctx) != -1 ||
            !fmpz_mpoly_q_factored_is_zero(A, ctx))
        {
            flint_printf("FAIL: division by zero\n");
            flint_abort();
        }

        if (fmpz_mpoly_q_factored_set_str_pretty(A, "x+z", vars, ctx) != -1)
        {
            flint_printf("FAIL: unknown variable\n");
            flint_abort();
        }

        fmpz_mpoly_q_factored_clear(A, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
