/*
    Copyright (C) 2010, 2012 William Hart
    Copyright (C) 2011, 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result = 1;
    fmpz_mod_ctx_t ctx;
    FLINT_TEST_INIT(state);

    flint_printf("evaluate_fmpz_vec_fast....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t P;
        fmpz * x, * y, * z;
        fmpz_t mod;
        slong j, n, npoints;

        fmpz_init(mod);
        
        do 
        {
           fmpz_randtest_unsigned(mod, state, 5);
           fmpz_add_ui(mod, mod, 2);
        } while (!fmpz_is_probabprime(mod));

        fmpz_mod_ctx_set_modulus(ctx, mod);

        npoints = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_mod_poly_init(P, ctx);
        x = _fmpz_vec_init(npoints);
        y = _fmpz_vec_init(npoints);
        z = _fmpz_vec_init(npoints);

        fmpz_mod_poly_randtest(P, state, n, ctx);

        for (j = 0; j < npoints; j++)
            fmpz_randtest_mod(x + j, state, mod);

        fmpz_mod_poly_evaluate_fmpz_vec_iter(y, P, x, npoints, ctx);
        fmpz_mod_poly_evaluate_fmpz_vec_fast(z, P, x, npoints, ctx);
        
        result = _fmpz_vec_equal(y, z, npoints);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("mod=");
            fmpz_print(mod);
            flint_printf(", n=%wd, npoints=%wd\n\n", n, npoints);
            flint_printf("P: "); fmpz_mod_poly_print(P, ctx); flint_printf("\n\n");
            for (j = 0; j < npoints; j++)
               fmpz_print(y + j), flint_printf(" ");
            flint_printf("\n");
            for (j = 0; j < npoints; j++)
               fmpz_print(z + j), flint_printf(" ");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(mod);
        fmpz_mod_poly_clear(P, ctx);
        _fmpz_vec_clear(x, npoints);
        _fmpz_vec_clear(y, npoints);
        _fmpz_vec_clear(z, npoints);
    }

    fmpz_mod_ctx_clear(ctx);
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
