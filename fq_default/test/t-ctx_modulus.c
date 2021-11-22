/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default.h"

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("ctx_modulus....");
    fflush(stdout);

    /* fq_zech range */
    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
	fmpz_mod_ctx_t mod_ctx;
        fmpz_mod_poly_t mod;
	fmpz_t p;

        fmpz_init(p);
        
        fmpz_set_ui(p, 3);

        fq_default_ctx_init(ctx, p, 3, "x");

	fmpz_mod_ctx_init(mod_ctx, p);
	fmpz_mod_poly_init(mod, mod_ctx);

        fq_default_ctx_modulus(mod, ctx);

        result = (fmpz_mod_poly_length(mod, mod_ctx) == 4);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(mod, mod_ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(mod, mod_ctx);
	fmpz_mod_ctx_clear(mod_ctx);

        fq_default_ctx_clear(ctx);
    }
    
    /* fq_nmod range */
    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fmpz_mod_ctx_t mod_ctx;
	fmpz_mod_poly_t mod;
        fmpz_t p;

        fmpz_init(p);

        fmpz_set_ui(p, 3);

        fq_default_ctx_init(ctx, p, 16, "x");

        fmpz_mod_ctx_init(mod_ctx, p);
        fmpz_mod_poly_init(mod, mod_ctx);

        fq_default_ctx_modulus(mod, ctx);

        result = (fmpz_mod_poly_length(mod, mod_ctx) == 17);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(mod, mod_ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(mod, mod_ctx);
        fmpz_mod_ctx_clear(mod_ctx);

        fq_default_ctx_clear(ctx);

        fmpz_clear(p);
    }

    /* fq range */
    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fmpz_mod_ctx_t mod_ctx;
	fmpz_mod_poly_t mod;
        fmpz_t p;

        fmpz_init(p);

        fmpz_set_str(p, "73786976294838206473", 10);

	fq_default_ctx_init(ctx, p, 1, "x");

        fmpz_mod_ctx_init(mod_ctx, p);
        fmpz_mod_poly_init(mod, mod_ctx);

        fq_default_ctx_modulus(mod, ctx);

        result = (fmpz_mod_poly_length(mod, mod_ctx) == 2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(mod, mod_ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(mod, mod_ctx);
        fmpz_mod_ctx_clear(mod_ctx);

        fq_default_ctx_clear(ctx);

        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
        return 0;
}
