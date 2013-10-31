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
    int i;
    flint_rand_t state;
    flint_randinit(state);
    printf("compose_mod_brent_kung_preinv....");
    fflush(stdout);

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, cinv, d, e;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, poly_init)(a, ctx);
        TEMPLATE(T, poly_init)(b, ctx);
        TEMPLATE(T, poly_init)(c, ctx);
        TEMPLATE(T, poly_init)(cinv, ctx);
        TEMPLATE(T, poly_init)(d, ctx);
        TEMPLATE(T, poly_init)(e, ctx);

        TEMPLATE(T, poly_randtest)(a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest)(b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero)(c, state, n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_reverse)(cinv, c, c->length, ctx);
        TEMPLATE(T, poly_inv_series_newton)(cinv, cinv, c->length, ctx);

        TEMPLATE(T, poly_rem)(a, a, c, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(d, a, b, c, cinv, ctx);
        TEMPLATE(T, poly_compose)(e, a, b, ctx);
        TEMPLATE(T, poly_rem)(e, e, c, ctx);

        if (!TEMPLATE(T, poly_equal)(d, e, ctx))
        {
            printf("FAIL (composition):\n");
            printf("a:\n"); TEMPLATE(T, poly_print)(a, ctx); printf("\n");
            printf("b:\n"); TEMPLATE(T, poly_print)(b, ctx); printf("\n");
            printf("c:\n"); TEMPLATE(T, poly_print)(c, ctx); printf("\n");
            printf("d:\n"); TEMPLATE(T, poly_print)(d, ctx); printf("\n");
            printf("e:\n"); TEMPLATE(T, poly_print)(e, ctx); printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, poly_clear)(b, ctx);
        TEMPLATE(T, poly_clear)(c, ctx);
        TEMPLATE(T, poly_clear)(cinv, ctx);
        TEMPLATE(T, poly_clear)(d, ctx);
        TEMPLATE(T, poly_clear)(e, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Test aliasing of res and a */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, cinv, d;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, poly_init)(a, ctx);
        TEMPLATE(T, poly_init)(b, ctx);
        TEMPLATE(T, poly_init)(c, ctx);
        TEMPLATE(T, poly_init)(cinv, ctx);
        TEMPLATE(T, poly_init)(d, ctx);

        TEMPLATE(T, poly_randtest)(a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest)(b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero)(c, state, n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_reverse)(cinv, c, c->length, ctx);
        TEMPLATE(T, poly_inv_series_newton)(cinv, cinv, c->length, ctx);

        TEMPLATE(T, poly_rem)(a, a, c, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(d, a, b, c, cinv, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(a, a, b, c, cinv, ctx);

        if (!TEMPLATE(T, poly_equal)(d, a, ctx))
        {
            printf("FAIL (aliasing a):\n");
            printf("a:\n"); TEMPLATE(T, poly_print)(a, ctx); printf("\n");
            printf("b:\n"); TEMPLATE(T, poly_print)(b, ctx); printf("\n");
            printf("c:\n"); TEMPLATE(T, poly_print)(c, ctx); printf("\n");
            printf("d:\n"); TEMPLATE(T, poly_print)(d, ctx); printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, poly_clear)(b, ctx);
        TEMPLATE(T, poly_clear)(c, ctx);
        TEMPLATE(T, poly_clear)(cinv, ctx);
        TEMPLATE(T, poly_clear)(d, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Test aliasing of res and b */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, cinv, d;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        
        TEMPLATE(T, poly_init)(a, ctx);
        TEMPLATE(T, poly_init)(b, ctx);
        TEMPLATE(T, poly_init)(c, ctx);
        TEMPLATE(T, poly_init)(cinv, ctx);
        TEMPLATE(T, poly_init)(d, ctx);

        TEMPLATE(T, poly_randtest)(a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest)(b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero)(c, state, n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_reverse)(cinv, c, c->length, ctx);
        TEMPLATE(T, poly_inv_series_newton)(cinv, cinv, c->length, ctx);

        TEMPLATE(T, poly_rem)(a, a, c, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(d, a, b, c, cinv, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(b, a, b, c, cinv, ctx);

        if (!TEMPLATE(T, poly_equal)(d, b, ctx))
        {
            printf("FAIL (aliasing b)\n");
            printf("a:\n"); TEMPLATE(T, poly_print)(a, ctx); printf("\n");
            printf("b:\n"); TEMPLATE(T, poly_print)(b, ctx); printf("\n");
            printf("c:\n"); TEMPLATE(T, poly_print)(c, ctx); printf("\n");
            printf("d:\n"); TEMPLATE(T, poly_print)(d, ctx); printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, poly_clear)(b, ctx);
        TEMPLATE(T, poly_clear)(c, ctx);
        TEMPLATE(T, poly_clear)(cinv, ctx);
        TEMPLATE(T, poly_clear)(d, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Test aliasing of res and c */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, cinv, d;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, poly_init)(a, ctx);
        TEMPLATE(T, poly_init)(b, ctx);
        TEMPLATE(T, poly_init)(c, ctx);
        TEMPLATE(T, poly_init)(cinv, ctx);
        TEMPLATE(T, poly_init)(d, ctx);

        TEMPLATE(T, poly_randtest)(a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest)(b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero)(c, state, n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_reverse)(cinv, c, c->length, ctx);
        TEMPLATE(T, poly_inv_series_newton)(cinv, cinv, c->length, ctx);

        TEMPLATE(T, poly_rem)(a, a, c, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(d, a, b, c, cinv, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(c, a, b, c, cinv, ctx);

        if (!TEMPLATE(T, poly_equal)(d, c, ctx))
        {
            printf("FAIL (aliasing c)\n");
            printf("a:\n"); TEMPLATE(T, poly_print)(a, ctx); printf("\n");
            printf("b:\n"); TEMPLATE(T, poly_print)(b, ctx); printf("\n");
            printf("c:\n"); TEMPLATE(T, poly_print)(c, ctx); printf("\n");
            printf("d:\n"); TEMPLATE(T, poly_print)(d, ctx); printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, poly_clear)(b, ctx);
        TEMPLATE(T, poly_clear)(c, ctx);
        TEMPLATE(T, poly_clear)(cinv, ctx);
        TEMPLATE(T, poly_clear)(d, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}


#endif
