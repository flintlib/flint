/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

#include <stdio.h>

#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("norm... ");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with product of Galois conjugates */
    for (i = 0; i < 2000; i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;
        fmpz_t x, y;
        long j;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);
        fmpz_init(x);
        fmpz_init(y);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, reduce)(a, ctx);

        TEMPLATE(T, norm)(x, a, ctx);

        TEMPLATE(T, one)(b, ctx);
        for (j = 0; j < TEMPLATE(T, ctx_degree)(ctx); j++)
        {
            TEMPLATE(T, frobenius)(c, a, j, ctx);
            TEMPLATE(T, mul)(b, b, c, ctx);
        }
        fmpz_poly_get_coeff_fmpz(y, b, 0);

        result = fmpz_equal(x, y);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("y = "), fmpz_print(y), flint_printf("\n");
            for (j = 0; j < TEMPLATE(T, ctx_degree)(ctx); j++)
            {
                TEMPLATE(T, frobenius)(c, a, j, ctx);
                flint_printf("sigma^%ld = ", j), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            }
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);
        fmpz_clear(x);
        fmpz_clear(y);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


#endif
