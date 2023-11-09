/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fmpz_mpoly_factor.h"

#define FLINT_ARRAY_ALLOC(n, T) (T *) flint_malloc((n)*sizeof(T))

/*
    If BB and CC are squarefree factorizations, then A should be as well.
*/
int fmpz_mpoly_factor_mul(
    fmpz_mpoly_factor_t A,
    const fmpz_mpoly_factor_t BB,
    const fmpz_mpoly_factor_t CC,
    const fmpz_mpoly_ctx_t ctx)
{
    int success = 0;
    slong i, j;
    fmpz_mpoly_factor_t B, C;
    fmpz_t t;

    fmpz_init(t);
    fmpz_mpoly_factor_init(B, ctx);
    fmpz_mpoly_factor_init(C, ctx);

    fmpz_mpoly_factor_set(B, BB, ctx);
    fmpz_mpoly_factor_set(C, CC, ctx);

    fmpz_mul(A->constant, B->constant, C->constant);
    A->num = 0;

    for (i = 0; i < B->num; i++)
    for (j = 0; j < C->num; j++)
    {
        fmpz_mpoly_factor_fit_length(A, A->num + 1, ctx);

        fmpz_add(A->exp + A->num, B->exp + i, C->exp + j);
        if (!fmpz_mpoly_gcd_cofactors(A->poly + A->num, B->poly + i,
                                   C->poly + j, B->poly + i, C->poly + j, ctx))
            goto cleanup;

        if (fmpz_is_zero(A->exp + A->num))
            continue;

        if (fmpz_mpoly_is_fmpz(A->poly + A->num, ctx))
        {
            fmpz_mpoly_get_fmpz(t, A->poly + A->num, ctx);
            if (!fmpz_pow_fmpz(t, t, A->exp + A->num))
                goto cleanup;
            fmpz_mul(A->constant, A->constant, t);
        }
        else
        {
            A->num++;
        }
    }

    fmpz_mpoly_factor_fit_length(A, A->num + B->num, ctx);
    for (i = 0; i < B->num; i++)
    {
        if (fmpz_is_zero(B->exp + i))
            continue;

        if (fmpz_mpoly_is_fmpz(B->poly + i, ctx))
        {
            fmpz_mpoly_get_fmpz(t, B->poly + i, ctx);
            if (!fmpz_pow_fmpz(t, t, B->exp + i))
                goto cleanup;
            fmpz_mul(A->constant, A->constant, t);
        }
        else
        {
            fmpz_mpoly_swap(A->poly + A->num, B->poly + i, ctx);
            fmpz_swap(A->exp + A->num, B->exp + i);
            A->num++;
        }
    }

    fmpz_mpoly_factor_fit_length(A, A->num + C->num, ctx);
    for (j = 0; j < C->num; j++)
    {
        if (fmpz_is_zero(C->exp + j))
            continue;

        if (fmpz_mpoly_is_fmpz(C->poly + j, ctx))
        {
            fmpz_mpoly_get_fmpz(t, C->poly + j, ctx);
            if (!fmpz_pow_fmpz(t, t, C->exp + j))
                goto cleanup;
            fmpz_mul(A->constant, A->constant, t);
        }
        else
        {
            fmpz_mpoly_swap(A->poly + A->num, C->poly + j, ctx);
            fmpz_swap(A->exp + A->num, C->exp + j);
            A->num++;
        }
    }

    success = 1;

cleanup:

    fmpz_clear(t);
    fmpz_mpoly_factor_clear(B, ctx);
    fmpz_mpoly_factor_clear(C, ctx);

    return success;
}

TEST_FUNCTION_START(fmpz_mpoly_factor_lcc_kaltofen, state)
{
    slong i, j, k, l;

    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        int had_zero;
        slong r, v, nvars;
        ulong c, * bounds;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_struct A[1], t1[1], t2[1], * divs, * lcs;
        fmpz_mpoly_factor_t Af, tf;
        fmpz_poly_struct ut1[1], ut2[1], * ulcs;
        fmpz * alphas, * content_divs;
        fmpz_t g1, g2, g3;

        nvars = 3 + n_randint(state, 4);

        fmpz_mpoly_ctx_init(ctx, nvars, ORD_LEX);
        fmpz_mpoly_init(A, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);
        fmpz_mpoly_factor_init(Af, ctx);
        fmpz_mpoly_factor_init(tf, ctx);
        fmpz_poly_init(ut1);
        fmpz_poly_init(ut2);

        fmpz_init(g1);
        fmpz_init(g2);
        fmpz_init(g3);

        bounds = FLINT_ARRAY_ALLOC(nvars, ulong);
        alphas = _fmpz_vec_init(nvars - 1);

        fmpz_mpoly_factor_one(Af, ctx);

        r = 2 + n_randint(state, 4);
        divs = FLINT_ARRAY_ALLOC(r, fmpz_mpoly_struct);
        lcs = FLINT_ARRAY_ALLOC(r, fmpz_mpoly_struct);
        ulcs = FLINT_ARRAY_ALLOC(r, fmpz_poly_struct);
        content_divs = _fmpz_vec_init(r);
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_init(divs + j, ctx);
            fmpz_mpoly_one(divs + j, ctx);
            fmpz_mpoly_init(lcs + j, ctx);
            c = 1 + n_randint(state, 100);
            fmpz_mul_ui(Af->constant, Af->constant, c);
            fmpz_mpoly_set_ui(lcs + j, c, ctx);
            fmpz_poly_init(ulcs + j);
        }

        for (j = 0; j < 1 + 2*r/3; j++)
        {
            bounds[0] = 1;
            for (k = 1; k < nvars; k++)
                bounds[k] = 1 + n_randint(state, 4);

            fmpz_mpoly_randtest_bounds(t1, state, 2 + 15/r, 20, bounds, ctx);
            if (fmpz_mpoly_is_zero(t1, ctx))
                fmpz_mpoly_one(t1, ctx);
            if (!fmpz_mpoly_factor_squarefree(tf, t1, ctx))
            {
                flint_printf("FAIL:\ncheck factor_squarefree success\n");
                fflush(stdout);
                flint_abort();
            }

            for (k = n_randint(state, 3); k >= 0; k--)
            {
                l = n_randint(state, r);
                fmpz_mpoly_mul(lcs + l, lcs + l, t1, ctx);
                fmpz_mpoly_factor_mul(Af, Af, tf, ctx);
            }
        }

        fmpz_mpoly_one(A, ctx);
        for (j = 0; j < r; j++)
            fmpz_mpoly_mul(A, A, lcs + j, ctx);

        for (j = 0; j < nvars - 1; j++)
        {
            fmpz_set_ui(alphas + j, n_urandint(state, 100));
            if (n_randint(state, 2))
                fmpz_neg(alphas + j, alphas + j);
        }

        had_zero = 0;
        for (v = 1; v < nvars; v++)
        {
            int have_zero = 0;

            for (j = 0; j < r; j++)
            {
                fmpz_mpoly_evaluate_rest_except_one(ut1, lcs + j, alphas, v, ctx);
                fmpz_mpoly_evaluate_rest_except_one(ut2, divs + j, alphas, v, ctx);
                if (fmpz_poly_is_zero(ut1) || fmpz_poly_is_zero(ut2))
                {
                    have_zero = 1;
                    had_zero = 1;
                }
                else
                {
                    fmpz_poly_primitive_part(ut2, ut2);
                    if (!fmpz_poly_divides(ulcs + j, ut1, ut2))
                    {
                        flint_printf("FAIL:\nbad divisor\n");
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }

            if (have_zero)
                continue;

            fmpz_mpoly_factor_lcc_kaltofen_step(divs, r, Af, ulcs, v, alphas, ctx);
        }

        if (Af->num == 0 && !had_zero)
        {
            for (j = 0; j < r; j++)
                fmpz_one(content_divs + j);

            for (v = 1; v < nvars; v++)
            {
                for (j = 0; j < r; j++)
                {
                    fmpz_mpoly_evaluate_rest_except_one(ut1, lcs + j, alphas, v, ctx);
                    fmpz_mpoly_evaluate_rest_except_one(ut2, divs + j, alphas, v, ctx);
                    _fmpz_vec_content(g1, ut1->coeffs, ut1->length);
                    _fmpz_vec_content(g2, ut2->coeffs, ut2->length);
                    fmpz_gcd(g3, g1, g2);
                    fmpz_divexact(g1, g1, g3);
                    fmpz_lcm(content_divs + j, content_divs + j, g1);
                }
            }

            for (j = 0; j < r; j++)
            {
                if (!fmpz_divisible(Af->constant, content_divs + j))
                {
                    flint_printf("FAIL:\nbad divisor\n");
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_divexact(Af->constant, Af->constant, content_divs + j);

                fmpz_mpoly_scalar_mul_fmpz(divs + j, divs + j, content_divs + j, ctx);

                if (!fmpz_mpoly_divides(t1, lcs + j, divs + j, ctx))
                {
                    flint_printf("FAIL:\nbad divisor\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        _fmpz_vec_clear(content_divs, r);
        fmpz_clear(g1);
        fmpz_clear(g2);
        fmpz_clear(g3);

        flint_free(bounds);
        _fmpz_vec_clear(alphas, nvars - 1);

        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_clear(divs + j, ctx);
            fmpz_mpoly_clear(lcs + j, ctx);
            fmpz_poly_clear(ulcs + j);
        }
        flint_free(divs);
        flint_free(lcs);
        flint_free(ulcs);

        fmpz_mpoly_clear(A, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
        fmpz_mpoly_factor_clear(Af, ctx);
        fmpz_mpoly_factor_clear(tf, ctx);
        fmpz_poly_clear(ut1);
        fmpz_poly_clear(ut2);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
