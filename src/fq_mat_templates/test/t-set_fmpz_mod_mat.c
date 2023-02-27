/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>

#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("set_fmpz_mod_mat... ");
    fflush(stdout);

    /* Check conversion of identity matrix */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) a;
        fmpz_mod_mat_t m;
        slong r, c;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        r = n_randint(state, 10);
	c = n_randint(state, 10);
	TEMPLATE(T, mat_init)(a, r, c, ctx);

        TEMPLATE(T, mat_randtest)(a, state, ctx);

        fmpz_mod_mat_init(m, r, c, TEMPLATE(T, ctx_prime)(ctx));

	fmpz_mod_mat_one(m);

        TEMPLATE(T, mat_set_fmpz_mod_mat)(a, m, ctx);

        result = (TEMPLATE(T, mat_is_one)(a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, mat_print)(a, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(m);

        TEMPLATE(T, mat_clear)(a, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}



#endif
