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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"


int
main(void)
{
    int iter;
    FLINT_TEST_INIT(state);

    flint_printf("randtest_irreducible....");
    fflush(stdout);

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) poly;
        slong length;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (poly, ctx);

        length = n_randint(state, 20) + 2;
        TEMPLATE(T, poly_randtest_irreducible) (poly, state, length, ctx);

        if (!TEMPLATE(T, poly_is_irreducible) (poly, ctx))
        {
            flint_printf("Error: reducible polynomial created!\n");
            flint_printf("poly:\n");
            TEMPLATE(T, poly_print_pretty) (poly, "x", ctx);
            flint_printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (poly, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}


#endif
