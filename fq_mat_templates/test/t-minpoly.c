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

    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2015 William Hart

******************************************************************************/


#ifdef T

#include "templates.h"

int
main(void)
{
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, mat_t) A;
    TEMPLATE(T, poly_t) p1, p2, q, r;
    slong i, m, n;
    FLINT_TEST_INIT(state);

    flint_printf("minpoly....");
    fflush(stdout);

    /* minpoly(A) divides charpoly(A) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);
        n = m;
        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (p1, ctx);
        TEMPLATE(T, poly_init) (p2, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_randtest) (A, state, ctx);
        
        TEMPLATE(T, mat_charpoly) (p1, A, ctx);
        TEMPLATE(T, mat_minpoly) (p2, A, ctx);
        
        TEMPLATE(T, poly_divrem) (q, r, p1, p2, ctx);
        
        if (!TEMPLATE(T, poly_is_zero) (r, ctx))
        {
            flint_printf("FAIL:\n");
            flint_printf("minpoly(A) doesn't divide charpoly(A).\n");
            abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, poly_clear) (p1, ctx);
        TEMPLATE(T, poly_clear) (p2, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}


#endif
