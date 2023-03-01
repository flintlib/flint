/*
    Copyright (C) 2011, 2021, 2022 William Hart
    Copyright (C) 2011 Fredrik Johansson

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
    
    flint_printf("sqrt_series....");
    fflush(stdout);

    /* Check g^2 = h mod x^m */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) h, g, r;
        slong m;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) one;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) != 0)
        {
            TEMPLATE(T, init)(one, ctx);

            TEMPLATE(T, poly_init)(h, ctx);
            TEMPLATE(T, poly_init)(g, ctx);
            TEMPLATE(T, poly_init)(r, ctx);

            do TEMPLATE(T, poly_randtest)(h, state, n_randint(state, 1000), ctx);
            while (h->length == 0);
            TEMPLATE(T, set_ui)(one, 1, ctx);
            TEMPLATE(T, poly_set_coeff)(h, 0, one, ctx);

            m = n_randint(state, h->length) + 1;

            TEMPLATE(T, poly_sqrt_series)(g, h, m, ctx);
            TEMPLATE(T, poly_mullow)(r, g, g, m, ctx);
            TEMPLATE(T, poly_truncate)(h, m, ctx);

            result = (TEMPLATE(T, poly_equal)(r, h, ctx));
            if (!result)
            {
                flint_printf("FAIL:\n");
                TEMPLATE(T, poly_print)(h, ctx), flint_printf("\n\n");
                TEMPLATE(T, poly_print)(g, ctx), flint_printf("\n\n");
                TEMPLATE(T, poly_print)(r, ctx), flint_printf("\n\n");
                flint_printf("n = ");
                fmpz_print(TEMPLATE(T, ctx_prime)(ctx));
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        
            TEMPLATE(T, poly_clear)(h, ctx);
            TEMPLATE(T, poly_clear)(g, ctx);
            TEMPLATE(T, poly_clear)(r, ctx);

            TEMPLATE(T, clear)(one, ctx);
       }

       TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check aliasing of h and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) g, h;
        slong m;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) one;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) != 0)
        {
            TEMPLATE(T, init)(one, ctx);

            TEMPLATE(T, poly_init)(h, ctx);
            TEMPLATE(T, poly_init)(g, ctx);
            
            do TEMPLATE(T, poly_randtest)(h, state, n_randint(state, 500), ctx);
            while (h->length == 0);
            TEMPLATE(T, set_ui)(one, 1, ctx);
            TEMPLATE(T, poly_set_coeff)(h, 0, one, ctx);

            m = n_randint(state, h->length) + 1;

            TEMPLATE(T, poly_sqrt_series)(g, h, m, ctx);
            TEMPLATE(T, poly_sqrt_series)(h, h, m, ctx);

            result = (TEMPLATE(T, poly_equal)(g, h, ctx));
            if (!result)
            {
                flint_printf("FAIL:\n");
                TEMPLATE(T, poly_print)(h, ctx), flint_printf("\n\n");
                TEMPLATE(T, poly_print)(g, ctx), flint_printf("\n\n");
                flint_printf("n = ");
                fmpz_print(TEMPLATE(T, ctx_prime)(ctx));
                flint_printf("m = %wd\n", m);
                fflush(stdout);
                flint_abort();
            }

            TEMPLATE(T, poly_clear)(g, ctx);
            TEMPLATE(T, poly_clear)(h, ctx);

            TEMPLATE(T, clear)(one, ctx);
        }

        TEMPLATE(T, ctx_clear)(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif

