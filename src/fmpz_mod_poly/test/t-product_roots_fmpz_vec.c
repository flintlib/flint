/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_product_roots_fmpz_vec, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t P, Q, tmp;
        fmpz * x;
        slong j, n, bits;
        fmpz_t mod;
        fmpz_mod_ctx_t ctx;

        n = n_randint(state, 100);
        bits = n_randint(state, 100);

        fmpz_init(mod);
        fmpz_randtest_unsigned(mod, state, bits);
        fmpz_add_ui(mod, mod, 2);
        fmpz_mod_ctx_init(ctx, mod);

        x = _fmpz_vec_init(n);
        _fmpz_vec_randtest(x, state, n, bits);

        for (j = 0; j < n; j++)
            fmpz_mod_set_fmpz(x + j, x + j, ctx);

        fmpz_mod_poly_init(P, ctx);
        fmpz_mod_poly_init(Q, ctx);
        fmpz_mod_poly_init(tmp, ctx);

        fmpz_mod_poly_product_roots_fmpz_vec(P, x, n, ctx);

        fmpz_mod_poly_one(Q, ctx);
        for (j = 0; j < n; j++)
        {
            fmpz_mod_poly_zero(tmp, ctx);
            fmpz_mod_poly_set_coeff_si(tmp, 1, -WORD(1), ctx);
            fmpz_mod_poly_set_coeff_fmpz(tmp, 0, x + j, ctx);
            fmpz_mod_poly_neg(tmp, tmp, ctx);
            fmpz_mod_poly_mul(Q, Q, tmp, ctx);
        }

        if (!fmpz_mod_poly_equal(P, Q, ctx))
        {
            flint_printf("FAIL (P != Q):\n");
            fmpz_mod_poly_print(P, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(Q, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(P, ctx);
        fmpz_mod_poly_clear(Q, ctx);
        fmpz_mod_poly_clear(tmp, ctx);
        _fmpz_vec_clear(x, n);
        fmpz_mod_ctx_clear(ctx);
        fmpz_clear(mod);
    }

    TEST_FUNCTION_END(state);
}
