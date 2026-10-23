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
#include "gr_vec.h"
#include "gr_poly.h"

/* Corner cases of modular powering (moduli of degree 0, 1 and 2, small
   exponents, aliasing, the wrappers taking a ulong exponent and the
   underscore functions) and of deflation. */
TEST_FUNCTION_START(gr_poly_powmod_misc, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t A, F, Finv, res1, res2;
        fmpz_t e;
        ulong ee;
        slong lenf;
        int status = GR_SUCCESS;

        gr_ctx_init_random_commutative_ring(ctx, state);

        gr_poly_init(A, ctx);
        gr_poly_init(F, ctx);
        gr_poly_init(Finv, ctx);
        gr_poly_init(res1, ctx);
        gr_poly_init(res2, ctx);
        fmpz_init(e);

        lenf = n_randint(state, 5);
        status |= gr_poly_randtest(F, state, lenf, ctx);
        if (lenf >= 1)
            status |= gr_poly_set_coeff_ui(F, lenf - 1, 1, ctx);
        status |= gr_poly_randtest(A, state, 1 + n_randint(state, 6), ctx);

        ee = n_randint(state, 5);
        fmpz_set_ui(e, ee);

        /* the fmpz and ulong versions agree */
        status |= gr_poly_powmod_fmpz_binexp(res1, A, e, F, ctx);
        status |= gr_poly_powmod_ui_binexp(res2, A, ee, F, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(res1, res2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (ui vs fmpz)\n\n");
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_printf("e = %wu\n", ee);
            flint_abort();
        }

        /* aliasing res = A */
        status |= gr_poly_set(res2, A, ctx);
        status |= gr_poly_powmod_ui_binexp(res2, res2, ee, F, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(res1, res2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (aliasing)\n\n");
            gr_ctx_println(ctx);
            flint_abort();
        }

        /* versions taking a precomputed inverse; for short moduli the
           inverse is not used, which is the case of interest here */
        if (lenf >= 1)
        {
            status |= gr_poly_reverse(Finv, F, lenf, ctx);
            status |= gr_poly_inv_series(Finv, Finv, lenf, ctx);

            if (status == GR_SUCCESS)
            {
                status |= gr_poly_powmod_ui_binexp_preinv(res2, A, ee, F, Finv, ctx);

                if (status == GR_SUCCESS && gr_poly_equal(res1, res2, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (preinv)\n\n");
                    gr_ctx_println(ctx);
                    flint_abort();
                }

                /* x^e mod F */
                status |= gr_poly_gen(res2, ctx);
                status |= gr_poly_powmod_fmpz_binexp(res1, res2, e, F, ctx);
                status |= gr_poly_powmod_x_fmpz_preinv(res2, e, F, Finv, ctx);

                if (status == GR_SUCCESS && gr_poly_equal(res1, res2, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (powmod_x)\n\n");
                    gr_ctx_println(ctx);
                    flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
                    flint_printf("e = %wu\n", ee);
                    flint_abort();
                }
            }
        }

        /* deflation with inflation 0 and 1 */
        status |= gr_poly_inflate(res1, A, 1, ctx);
        if (status == GR_SUCCESS && gr_poly_equal(res1, A, ctx) == T_FALSE)
        {
            flint_printf("FAIL (inflate by 1)\n\n");
            gr_ctx_println(ctx);
            flint_abort();
        }

        status |= gr_poly_inflate(res1, A, 0, ctx);
        {
            gr_ptr t;
            GR_TMP_INIT(t, ctx);
            status |= gr_one(t, ctx);
            status |= gr_poly_evaluate(t, A, t, ctx);
            status |= gr_poly_set_scalar(res2, t, ctx);
            GR_TMP_CLEAR(t, ctx);
        }

        if (status == GR_SUCCESS && gr_poly_equal(res1, res2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (inflate by 0)\n\n");
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(F, ctx);
        gr_poly_clear(Finv, ctx);
        gr_poly_clear(res1, ctx);
        gr_poly_clear(res2, ctx);
        fmpz_clear(e);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
