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

TEST_FUNCTION_START(gr_poly_preinv, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t A, B, F, Q, R, Q2, R2, T;
        gr_poly_preinv_t P;
        fmpz_t e;
        slong i, lenf;
        int which;
        int status = GR_SUCCESS;

        int finite = n_randint(state, 4) != 0;

        if (finite)
            gr_ctx_init_random_finite_field(ctx, state);
        else
            gr_ctx_init_fmpq(ctx);   /* exact, non-finite: coefficient growth */

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(F, ctx);
        gr_poly_init(Q, ctx);
        gr_poly_init(R, ctx);
        gr_poly_init(Q2, ctx);
        gr_poly_init(R2, ctx);
        gr_poly_init(T, ctx);
        gr_poly_preinv_init(P, ctx);
        fmpz_init(e);

        /* modulus: random, or sparse (few random terms), or with all
           coefficients 1 (the special case in the nmod loop) */
        lenf = 2 + n_randint(state, (finite && n_randint(state, 4) == 0) ? 100 : 12);
        which = n_randint(state, 3);

        status |= gr_poly_zero(F, ctx);
        if (which == 0)
        {
            status |= gr_poly_randtest(F, state, lenf - 1, ctx);
        }
        else
        {
            slong nterms = 1 + n_randint(state, 5);
            for (i = 0; i < nterms; i++)
            {
                if (which == 1)
                    status |= gr_poly_set_coeff_ui(F, n_randint(state, lenf - 1), 1, ctx);
                else
                {
                    gr_ptr c;
                    GR_TMP_INIT(c, ctx);
                    status |= gr_randtest(c, state, ctx);
                    status |= gr_poly_set_coeff_scalar(F, n_randint(state, lenf - 1), c, ctx);
                    GR_TMP_CLEAR(c, ctx);
                }
            }
        }

        /* leading coefficient: 1 or random invertible */
        if (n_randint(state, 2))
            status |= gr_poly_set_coeff_ui(F, lenf - 1, 1, ctx);
        else
        {
            gr_ptr c;
            GR_TMP_INIT(c, ctx);
            do {
                status |= gr_randtest(c, state, ctx);
            } while (gr_is_zero(c, ctx) != T_FALSE);
            status |= gr_poly_set_coeff_scalar(F, lenf - 1, c, ctx);
            GR_TMP_CLEAR(c, ctx);
        }

        status |= gr_poly_randtest(A, state, n_randint(state, 3 * lenf), ctx);
        status |= gr_poly_randtest(B, state, n_randint(state, lenf), ctx);

        switch (n_randint(state, 4))
        {
            case 0: status |= gr_poly_preinv_set(P, F, ctx); break;
            case 1: status |= gr_poly_preinv_set_newton(P, F, ctx); break;
            case 2: status |= gr_poly_preinv_set_sparse(P, F, ctx); break;
            default:
                /* transformed representation: not available for all
                   rings (and sizes); fall back to Newton */
                status |= gr_poly_preinv_set_transformed(P, F, ctx);
                if (status == GR_UNABLE)
                    status = gr_poly_preinv_set_newton(P, F, ctx);
                break;
        }

        if (status != GR_SUCCESS)
        {
            /* can fail over inexact or non-field rings; skip */
            goto cleanup;
        }

        /* divrem */
        status = gr_poly_preinv_divrem(Q, R, A, P, ctx);
        status |= gr_poly_divrem(Q2, R2, A, F, ctx);

        if (status == GR_SUCCESS && (gr_poly_equal(Q, Q2, ctx) == T_FALSE || gr_poly_equal(R, R2, ctx) == T_FALSE))
        {
            flint_printf("FAIL (divrem)\n\n");
            gr_ctx_println(ctx);
            flint_printf("kind = %d\n", P->kind);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_printf("Q = "); gr_poly_print(Q, ctx); flint_printf("\n");
            flint_printf("Q2 = "); gr_poly_print(Q2, ctx); flint_printf("\n");
            flint_printf("R = "); gr_poly_print(R, ctx); flint_printf("\n");
            flint_printf("R2 = "); gr_poly_print(R2, ctx); flint_printf("\n");
            flint_abort();
        }

        /* rem, aliased */
        status = gr_poly_set(R, A, ctx);
        status |= gr_poly_preinv_rem(R, R, P, ctx);
        if (status == GR_SUCCESS && gr_poly_equal(R, R2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (rem)\n\n");
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_abort();
        }

        /* mulmod */
        status = gr_poly_preinv_mulmod(T, A, B, P, ctx);
        status |= gr_poly_mulmod(R2, A, B, F, ctx);
        if (status == GR_SUCCESS && gr_poly_equal(T, R2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (mulmod)\n\n");
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_abort();
        }

        /* powmod: binexp, sliding, x */
        fmpz_randtest_unsigned(e, state, finite ? 60 : 5);
        status = gr_poly_preinv_powmod_fmpz_binexp(T, B, e, P, ctx);
        status |= gr_poly_powmod_fmpz_binexp(R2, B, e, F, ctx);
        if (status == GR_SUCCESS && gr_poly_equal(T, R2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (powmod binexp)\n\n");
            gr_ctx_println(ctx);
            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_printf("e = %{fmpz}\n", e);
            flint_abort();
        }

        status = gr_poly_preinv_powmod_fmpz_sliding(T, B, e, n_randint(state, 4), P, ctx);
        if (status == GR_SUCCESS && gr_poly_equal(T, R2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (powmod sliding)\n\n");
            gr_ctx_println(ctx);
            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_printf("e = %{fmpz}\n", e);
            flint_abort();
        }

        status = gr_poly_preinv_powmod_x_fmpz(T, e, P, ctx);
        status |= gr_poly_gen(B, ctx);
        status |= gr_poly_powmod_fmpz_binexp(R2, B, e, F, ctx);
        if (status == GR_SUCCESS && gr_poly_equal(T, R2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (powmod x)\n\n");
            gr_ctx_println(ctx);
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_printf("e = %{fmpz}\n", e);
            flint_printf("T = "); gr_poly_print(T, ctx); flint_printf("\n");
            flint_printf("R2 = "); gr_poly_print(R2, ctx); flint_printf("\n");
            flint_abort();
        }

cleanup:
        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(F, ctx);
        gr_poly_clear(Q, ctx);
        gr_poly_clear(R, ctx);
        gr_poly_clear(Q2, ctx);
        gr_poly_clear(R2, ctx);
        gr_poly_clear(T, ctx);
        gr_poly_preinv_clear(P, ctx);
        fmpz_clear(e);
        gr_ctx_clear(ctx);
    }

    /* Long dense monic moduli over nmod and mpn_mod: exercises the
       transformed representation (fft_small) when available. */
    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t A, B, F, Q, R, Q2, R2, T;
        gr_poly_preinv_t P;
        slong lenf;
        int status = GR_SUCCESS;

        if (n_randint(state, 2))
        {
            gr_ctx_init_nmod(ctx, n_randtest_prime(state, 1));
        }
        else
        {
            fmpz_t m;
            fmpz_init(m);
            fmpz_randprime(m, state, FLINT_BITS + 1 + n_randint(state, 2 * FLINT_BITS), 0);
            GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, m));
            fmpz_clear(m);
        }
        gr_poly_init(A, ctx); gr_poly_init(B, ctx); gr_poly_init(F, ctx);
        gr_poly_init(Q, ctx); gr_poly_init(R, ctx); gr_poly_init(Q2, ctx); gr_poly_init(R2, ctx); gr_poly_init(T, ctx);
        gr_poly_preinv_init(P, ctx);

        lenf = 500 + n_randint(state, 900);
        status |= gr_poly_randtest(F, state, lenf - 1, ctx);
        status |= gr_poly_set_coeff_ui(F, lenf - 1, 1, ctx);
        status |= gr_poly_randtest(A, state, n_randint(state, 3 * lenf), ctx);
        status |= gr_poly_randtest(B, state, n_randint(state, lenf), ctx);
        status |= gr_poly_preinv_set(P, F, ctx);

        status |= gr_poly_preinv_divrem(Q, R, A, P, ctx);
        status |= gr_poly_divrem(Q2, R2, A, F, ctx);
        status |= gr_poly_preinv_mulmod(T, A, B, P, ctx);

        if (status != GR_SUCCESS || gr_poly_equal(Q, Q2, ctx) == T_FALSE || gr_poly_equal(R, R2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (long modulus, divrem)\n\n");
            gr_ctx_println(ctx);
            flint_printf("kind = %d, lenf = %wd, lenA = %wd\n", P->kind, lenf, A->length);
            flint_abort();
        }

        status |= gr_poly_mulmod(R2, A, B, F, ctx);
        if (status != GR_SUCCESS || gr_poly_equal(T, R2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (long modulus, mulmod)\n\n");
            gr_ctx_println(ctx);
            flint_printf("kind = %d, lenf = %wd\n", P->kind, lenf);
            flint_abort();
        }

        gr_poly_clear(A, ctx); gr_poly_clear(B, ctx); gr_poly_clear(F, ctx);
        gr_poly_clear(Q, ctx); gr_poly_clear(R, ctx); gr_poly_clear(Q2, ctx); gr_poly_clear(R2, ctx); gr_poly_clear(T, ctx);
        gr_poly_preinv_clear(P, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
