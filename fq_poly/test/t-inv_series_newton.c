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

    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include "fq_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("inv_series_newton....");
    fflush(stdout);

    flint_randinit(state);

    /* Check Q^{-1} * Q is congruent 1 mod t^n */
    for (i = 0; i < 1000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, c, one;
        slong n = n_randint(state, 80) + 1;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(one, ctx);

        fq_poly_randtest_not_zero(a, state, n_randint(state, 80) + 1, ctx);
        fq_randtest_not_zero(a->coeffs, state, ctx);
        fq_poly_one(one, ctx);

        fq_poly_inv_series_newton(b, a, n, ctx);
        fq_poly_mullow(c, a, b, n, ctx);

        result = (fq_poly_equal(c, one, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fq_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c = "), fq_poly_print(c, ctx), flint_printf("\n\n");
            flint_printf("ctx = "), fq_ctx_print(ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(one, ctx);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

