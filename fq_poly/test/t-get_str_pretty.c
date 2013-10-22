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

#include "fq_poly.h"

int
main(void)
{
    int i, len;
    char *str;
    fq_poly_t a;
    fq_ctx_t ctx;
    flint_rand_t state;

    flint_printf("get_str_pretty....");
    fflush(stdout);

    flint_randinit(state);

    fq_ctx_randtest(ctx, state);
    
    fq_poly_init(a, ctx);
    for (len = 0; len < 100; len++)
        for (i = 0; i < 10; i++)
        {
            fq_poly_randtest(a, state, len, ctx);
            str = fq_poly_get_str_pretty(a, "x", ctx);
            /* printf("\n\n"); */
            /* fq_poly_print_pretty(a, "x", ctx); */
            /* printf("\n%s\n", str); */
            flint_free(str);
        }

    fq_poly_clear(a, ctx);
    fq_ctx_clear(ctx);

    flint_printf("PASS\n");
    return 0;
}
