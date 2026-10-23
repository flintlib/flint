/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(gr_poly_is_irreducible, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t A, B, P;
        gr_poly_vec_t fac;
        fmpz_vec_t exp;
        gr_ptr c;
        truth_t irr, irr_ddf, irr_ben_or, sqf, expected_irr, expected_sqf;
        slong i;
        int status = GR_SUCCESS;

        gr_ctx_init_random_finite_field(ctx, state);

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(P, ctx);
        gr_poly_vec_init(fac, 0, ctx);
        fmpz_vec_init(exp, 0);
        c = gr_heap_init(ctx);

        switch (n_randint(state, 3))
        {
            case 0:
                status |= gr_poly_randtest(P, state, 1 + n_randint(state, 30), ctx);
                break;
            case 1:
                /* product of two polynomials */
                status |= gr_poly_randtest(A, state, 1 + n_randint(state, 10), ctx);
                status |= gr_poly_randtest(B, state, 1 + n_randint(state, 10), ctx);
                status |= gr_poly_mul(P, A, B, ctx);
                break;
            default:
                /* square times something */
                status |= gr_poly_randtest(A, state, 1 + n_randint(state, 8), ctx);
                status |= gr_poly_randtest(B, state, 1 + n_randint(state, 8), ctx);
                status |= gr_poly_mul(P, A, A, ctx);
                status |= gr_poly_mul(P, P, B, ctx);
                break;
        }

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL (setup)\n");
            gr_ctx_println(ctx);
            flint_abort();
        }

        irr = gr_poly_is_irreducible(P, ctx);
        irr_ddf = gr_poly_is_irreducible_ddf(P, ctx);
        irr_ben_or = gr_poly_is_irreducible_ben_or(P, ctx);
        sqf = gr_poly_is_squarefree(P, ctx);

        if (P->length == 0)
        {
            expected_irr = T_FALSE;
            expected_sqf = T_FALSE;
        }
        else
        {
            status = gr_poly_factor_finite_field(c, fac, exp, P, ctx);

            if (status != GR_SUCCESS)
            {
                flint_printf("FAIL (factor)\n");
                gr_ctx_println(ctx);
                flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                flint_abort();
            }

            expected_irr = (fac->length == 1 && fmpz_is_one(exp->entries)) ? T_TRUE : T_FALSE;

            expected_sqf = T_TRUE;
            for (i = 0; i < fac->length; i++)
                if (!fmpz_is_one(exp->entries + i))
                    expected_sqf = T_FALSE;
        }

        if (irr != expected_irr || irr_ddf != expected_irr || irr_ben_or != expected_irr)
        {
            flint_printf("FAIL (is_irreducible)\n");
            gr_ctx_println(ctx);
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_printf("irr = %d, ddf = %d, ben_or = %d, expected = %d\n", irr, irr_ddf, irr_ben_or, expected_irr);
            flint_printf("fac = %{gr_poly*}\n", fac->entries, fac->length, ctx);
            flint_printf("exp = %{fmpz*}\n", exp->entries, exp->length);
            flint_abort();
        }

        if (sqf != expected_sqf)
        {
            flint_printf("FAIL (is_squarefree)\n");
            gr_ctx_println(ctx);
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_printf("sqf = %d, expected = %d\n", sqf, expected_sqf);
            flint_printf("fac = %{gr_poly*}\n", fac->entries, fac->length, ctx);
            flint_printf("exp = %{fmpz*}\n", exp->entries, exp->length);
            flint_abort();
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(P, ctx);
        gr_poly_vec_clear(fac, ctx);
        fmpz_vec_clear(exp);
        gr_heap_clear(c, ctx);
        gr_ctx_clear(ctx);
    }

    /* Large degrees, to reach the code path using Rabin's test:
       irreducible polynomials from nmod_poly_minimal_irreducible, and
       products of two of them. Two cases suffice to cover the path
       (irreducible and reducible input); the cost grows quickly with
       the degree, so the test multiplier only varies the field. */
    for (iter = 0; iter < 2; iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t P, Q;
        nmod_poly_t f;
        slong n, m, i;
        ulong p;
        truth_t r1, r2, r3;
        int reducible = iter;

        p = n_randprime(state, 2 + n_randint(state, 20), 1);
        gr_ctx_init_nmod(ctx, p);
        GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));
        gr_poly_init(P, ctx);
        gr_poly_init(Q, ctx);
        nmod_poly_init(f, p);

        n = 601 + n_randint(state, 20);   /* above the cutoff for Rabin's test */
        nmod_poly_minimal_irreducible(f, n);
        for (i = 0; i < f->length; i++)
            GR_MUST_SUCCEED(gr_poly_set_coeff_ui(P, i, f->coeffs[i], ctx));

        if (reducible)
        {
            m = 5 + n_randint(state, 20);
            nmod_poly_minimal_irreducible(f, m);
            for (i = 0; i < f->length; i++)
                GR_MUST_SUCCEED(gr_poly_set_coeff_ui(Q, i, f->coeffs[i], ctx));
            GR_MUST_SUCCEED(gr_poly_mul(P, P, Q, ctx));
        }

        /* dense: shift by 1 (preserves irreducibility) */
        if (n_randint(state, 2))
        {
            gr_ptr c;
            GR_TMP_INIT(c, ctx);
            GR_MUST_SUCCEED(gr_one(c, ctx));
            GR_MUST_SUCCEED(gr_poly_taylor_shift(P, P, c, ctx));
            GR_TMP_CLEAR(c, ctx);
        }

        r1 = gr_poly_is_irreducible(P, ctx);
        r2 = gr_poly_is_irreducible_rabin(P, ctx);
        r3 = gr_poly_is_irreducible_ddf(P, ctx);

        if (r1 != (reducible ? T_FALSE : T_TRUE) || r2 != r1 || r3 != r1)
        {
            flint_printf("FAIL (large degree)\n\n");
            gr_ctx_println(ctx);
            flint_printf("n = %wd, reducible = %d: irr = %d, rabin = %d, ddf = %d\n", n, reducible, r1, r2, r3);
            flint_abort();
        }

        gr_poly_clear(P, ctx);
        gr_poly_clear(Q, ctx);
        nmod_poly_clear(f);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
