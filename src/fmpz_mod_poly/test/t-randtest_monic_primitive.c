/*
    Copyright (C) 2018 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_randtest_monic_primitive, state)
{
    int i;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t poly;
        slong d;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state,
                                   FLINT_MIN(FLINT_BITS - 1, 10)), 1));
        fmpz_mod_ctx_set_modulus(ctx, p);
        d = n_randint(state, 5) + 3;

        fmpz_mod_poly_init(poly, ctx);
        fmpz_mod_poly_randtest_monic_primitive(poly, state, d, ctx);

        fmpz_mod_poly_clear(poly, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
