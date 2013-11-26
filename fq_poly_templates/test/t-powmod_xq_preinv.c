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
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("powmod_xq_preinv...");
    fflush(stdout);

    /* No aliasing */
    for (i = 0; i < 50; i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, res1, res2, t, f, finv;
        fmpz_t q;

        TEMPLATE(T, ctx_randtest)(ctx, state);


        TEMPLATE(T, poly_init)(a, ctx);
        TEMPLATE(T, poly_init)(f, ctx);
        TEMPLATE(T, poly_init)(finv, ctx);
        TEMPLATE(T, poly_init)(res1, ctx);
        TEMPLATE(T, poly_init)(res2, ctx);
        TEMPLATE(T, poly_init)(t, ctx);

        TEMPLATE(T, poly_gen)(a, ctx);
        TEMPLATE(T, poly_randtest_not_zero)(f, state, n_randint(state, 50) + 1, ctx);

        fmpz_init(q);
        TEMPLATE(T, ctx_order)(q, ctx);

        TEMPLATE(T, poly_reverse)(finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton)(finv, finv, f->length, ctx);

        TEMPLATE(T, poly_powmod_fmpz_binexp_preinv)(res1, a, q, f, finv, ctx);
        TEMPLATE(T, poly_powmod_xq_preinv)(res2, f, finv, ctx);

        result = (TEMPLATE(T, poly_equal)(res1, res2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, ctx_print)(ctx);
            flint_printf("f:\n"); TEMPLATE(T, poly_print_pretty)(f, "x", ctx);
            flint_printf("\n\n");
            flint_printf("res1:\n"); TEMPLATE(T, poly_print_pretty)(res1, "x", ctx);
            flint_printf("\n\n");
            flint_printf("res2:\n"); TEMPLATE(T, poly_print_pretty)(res2, "x", ctx);
            flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, poly_clear)(f, ctx);
        TEMPLATE(T, poly_clear)(finv, ctx);
        TEMPLATE(T, poly_clear)(res1, ctx);
        TEMPLATE(T, poly_clear)(res2, ctx);
        TEMPLATE(T, poly_clear)(t, ctx);
        fmpz_clear(q);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}


#endif
