/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_mulmid, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Compare with mullow */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;
        slong nlo, nhi;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        nlo = n_randint(state, 50);
        nhi = n_randint(state, 50);
        fmpz_mod_poly_randtest(b, state, 1 + n_randint(state, 50), ctx);
        fmpz_mod_poly_randtest(c, state, 1 + n_randint(state, 50), ctx);

        if (n_randint(state, 2))
        {
            fmpz_mod_poly_set(a, b, ctx);
            fmpz_mod_poly_mulmid(a, a, c, nlo, nhi, ctx);
        }
        else if (n_randint(state, 2))
        {
            fmpz_mod_poly_set(a, c, ctx);
            fmpz_mod_poly_mulmid(a, b, a, nlo, nhi, ctx);
        }
        else
        {
            fmpz_mod_poly_mulmid(a, b, c, nlo, nhi, ctx);
        }

        fmpz_mod_poly_mullow(b, b, c, nhi, ctx);
        fmpz_mod_poly_shift_right(b, b, nlo, ctx);

        result = (fmpz_mod_poly_equal(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
