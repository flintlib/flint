/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
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

    flint_printf("set_fmpz_poly... ");
    fflush(stdout);

    /* Check litfed polynomials by evaluating at random points */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) r, s;
        TEMPLATE(T, poly_t) a;
        fmpz_mod_poly_t b;
        fmpz_t p;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest)(ctx, state);
        TEMPLATE(T, init)(r, ctx); TEMPLATE(T, init)(s, ctx);
        TEMPLATE(T, poly_init)(a, ctx);
        fmpz_mod_poly_init(b, &(ctx->p));
        fmpz_init(p);

        fmpz_mod_poly_randtest(b, state, len);
        fmpz_randtest(p, state, 10);
        
        TEMPLATE(T, poly_set_fmpz_mod_poly)(a, b, ctx);
        TEMPLATE(T, set_fmpz)(r, p, ctx);
        TEMPLATE3(T, poly_evaluate, T)(r, a, r, ctx);
        fmpz_mod_poly_evaluate_fmpz(p, b, p);
        TEMPLATE(T, set_fmpz)(s, p, ctx);

        result = TEMPLATE(T, equal)(r, s, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), TEMPLATE(T, ctx_print(ctx)),
                flint_printf("\n");
            flint_printf("b = "), fmpz_mod_poly_print_pretty(b, "X"),
                flint_printf("\n");
            flint_printf("p = "), fmpz_print(p),
                flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(r, ctx); TEMPLATE(T, clear)(s, ctx);
        fmpz_mod_poly_clear(b);
        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

#endif
