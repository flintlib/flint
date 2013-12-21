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
    int i, len;
    char *str;
    TEMPLATE(T, poly_t) a;
    TEMPLATE(T, ctx_t) ctx;
    FLINT_TEST_INIT(state);

    flint_printf("get_str....");
    fflush(stdout);

    TEMPLATE(T, ctx_randtest) (ctx, state);

    TEMPLATE(T, poly_init) (a, ctx);
    for (len = 0; len < 100; len++)
        for (i = 0; i < 10; i++)
        {
            TEMPLATE(T, poly_randtest) (a, state, len, ctx);
            str = TEMPLATE(T, poly_get_str) (a, ctx);
            /* flint_printf("\n\n"); */
            /* TEMPLATE(T, poly_print)(a, ctx); */
            /* flint_printf("\n%s\n", str); */
            flint_free(str);
        }

    TEMPLATE(T, poly_clear) (a, ctx);
    TEMPLATE(T, ctx_clear) (ctx);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}


#endif
