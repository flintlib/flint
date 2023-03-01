/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    fmpz_mod_ctx_t ctx;
    FLINT_TEST_INIT(state);

    flint_printf("get/set_nmod_poly....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b;
        nmod_poly_t c;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, FLINT_BITS - 1);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        
        nmod_poly_init(c, fmpz_get_ui(p));

        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_get_nmod_poly(c, a);
        fmpz_mod_poly_set_nmod_poly(b, c);

        result = fmpz_mod_poly_equal(a, b, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n");
            flint_printf("c = "), nmod_poly_print(c), flint_printf("\n");
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        nmod_poly_clear(c);

        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

