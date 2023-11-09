/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

static int _check(fmpz_mod_poly_struct **B,
                  const fmpz_mod_poly_t F, const fmpz_mod_poly_t R,
                                            const fmpz_mod_ctx_t ctx)
{
    const slong lenF = F->length;
    const slong lenR = R->length;
    const slong N = (lenF - 1) / (lenR - 1);

    slong i;
    int result;

    fmpz_mod_poly_t S;

    fmpz_mod_poly_init(S, ctx);
    fmpz_mod_poly_set(S, B[N], ctx);
    for (i = N; i > 0; i--)
    {
        fmpz_mod_poly_mul(S, S, R, ctx);
        fmpz_mod_poly_add(S, S, B[i - 1], ctx);
    }
    result = fmpz_mod_poly_equal(F, S, ctx);
    fmpz_mod_poly_clear(S, ctx);
    return result;
}

TEST_FUNCTION_START(fmpz_mod_poly_radix, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t f, r;
        fmpz_mod_poly_struct **b;
        fmpz_mod_poly_radix_t D;
        slong j, N;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(r, ctx);
        fmpz_mod_poly_randtest(f, state, n_randint(state, 500), ctx);
        do
            fmpz_mod_poly_randtest_not_zero(r, state, n_randint(state, 20) + 2, ctx);
        while (r->length < 2);

        N = FLINT_MAX(0, fmpz_mod_poly_degree(f, ctx) / fmpz_mod_poly_degree(r, ctx));
        b = flint_malloc((N + 1) * sizeof(fmpz_mod_poly_struct *));
        for (j = 0; j <= N; j++)
        {
            b[j] = flint_malloc(sizeof(fmpz_mod_poly_struct));
            fmpz_mod_poly_init(b[j], ctx);
        }

        /* Ensure that lead(r) is a unit mod p */
        {
            fmpz_t d;
            fmpz *leadR = fmpz_mod_poly_lead(r, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadR);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadR, leadR, d);
                fmpz_gcd(d, p, leadR);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_radix_init(D, r, f->length - 1 + n_randint(state, 50), ctx);

        fmpz_mod_poly_radix(b, f, D, ctx);

        result = _check(b, f, r, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("result = %d\n", result);
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("f = "), fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            flint_printf("N = %wd\n\n", N);
            for (j = 0; j <= N; j++)
            {
                flint_printf("b[%wd] = ", j), fmpz_mod_poly_print(b[j], ctx), flint_printf("\n\n");
            }
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_mod_poly_radix_clear(D);
        for (j = 0; j <= N; j++)
        {
            fmpz_mod_poly_clear(b[j], ctx);
            flint_free(b[j]);
        }
        flint_free(b);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
