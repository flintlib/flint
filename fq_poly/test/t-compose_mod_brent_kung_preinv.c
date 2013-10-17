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
    int i;
    flint_rand_t state;
    flint_randinit(state);
    printf("compose_mod_brent_kung_preinv....");
    fflush(stdout);

    for (i = 0; i < 1000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, c, cinv, d, e;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(cinv, ctx);
        fq_poly_init(d, ctx);
        fq_poly_init(e, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fq_poly_reverse(cinv, c, c->length, ctx);
        fq_poly_inv_series_newton(cinv, cinv, c->length, ctx);

        fq_poly_rem(a, a, c, ctx);
        fq_poly_compose_mod_brent_kung_preinv(d, a, b, c, cinv, ctx);
        fq_poly_compose(e, a, b, ctx);
        fq_poly_rem(e, e, c, ctx);

        if (!fq_poly_equal(d, e, ctx))
        {
            printf("FAIL (composition):\n");
            printf("a:\n"); fq_poly_print(a, ctx); printf("\n");
            printf("b:\n"); fq_poly_print(b, ctx); printf("\n");
            printf("c:\n"); fq_poly_print(c, ctx); printf("\n");
            printf("d:\n"); fq_poly_print(d, ctx); printf("\n");
            printf("e:\n"); fq_poly_print(e, ctx); printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(cinv, ctx);
        fq_poly_clear(d, ctx);
        fq_poly_clear(e, ctx);

        fq_ctx_clear(ctx);
    }

    /* Test aliasing of res and a */
    for (i = 0; i < 1000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, c, cinv, d;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(cinv, ctx);
        fq_poly_init(d, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fq_poly_reverse(cinv, c, c->length, ctx);
        fq_poly_inv_series_newton(cinv, cinv, c->length, ctx);

        fq_poly_rem(a, a, c, ctx);
        fq_poly_compose_mod_brent_kung_preinv(d, a, b, c, cinv, ctx);
        fq_poly_compose_mod_brent_kung_preinv(a, a, b, c, cinv, ctx);

        if (!fq_poly_equal(d, a, ctx))
        {
            printf("FAIL (aliasing a):\n");
            printf("a:\n"); fq_poly_print(a, ctx); printf("\n");
            printf("b:\n"); fq_poly_print(b, ctx); printf("\n");
            printf("c:\n"); fq_poly_print(c, ctx); printf("\n");
            printf("d:\n"); fq_poly_print(d, ctx); printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(cinv, ctx);
        fq_poly_clear(d, ctx);

        fq_ctx_clear(ctx);
    }

    /* Test aliasing of res and b */
    for (i = 0; i < 1000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, c, cinv, d;

        fq_ctx_randtest(ctx, state);
        
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(cinv, ctx);
        fq_poly_init(d, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fq_poly_reverse(cinv, c, c->length, ctx);
        fq_poly_inv_series_newton(cinv, cinv, c->length, ctx);

        fq_poly_rem(a, a, c, ctx);
        fq_poly_compose_mod_brent_kung_preinv(d, a, b, c, cinv, ctx);
        fq_poly_compose_mod_brent_kung_preinv(b, a, b, c, cinv, ctx);

        if (!fq_poly_equal(d, b, ctx))
        {
            printf("FAIL (aliasing b)\n");
            printf("a:\n"); fq_poly_print(a, ctx); printf("\n");
            printf("b:\n"); fq_poly_print(b, ctx); printf("\n");
            printf("c:\n"); fq_poly_print(c, ctx); printf("\n");
            printf("d:\n"); fq_poly_print(d, ctx); printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(cinv, ctx);
        fq_poly_clear(d, ctx);

        fq_ctx_clear(ctx);
    }

    /* Test aliasing of res and c */
    for (i = 0; i < 1000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, c, cinv, d;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(cinv, ctx);
        fq_poly_init(d, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fq_poly_reverse(cinv, c, c->length, ctx);
        fq_poly_inv_series_newton(cinv, cinv, c->length, ctx);

        fq_poly_rem(a, a, c, ctx);
        fq_poly_compose_mod_brent_kung_preinv(d, a, b, c, cinv, ctx);
        fq_poly_compose_mod_brent_kung_preinv(c, a, b, c, cinv, ctx);

        if (!fq_poly_equal(d, c, ctx))
        {
            printf("FAIL (aliasing c)\n");
            printf("a:\n"); fq_poly_print(a, ctx); printf("\n");
            printf("b:\n"); fq_poly_print(b, ctx); printf("\n");
            printf("c:\n"); fq_poly_print(c, ctx); printf("\n");
            printf("d:\n"); fq_poly_print(d, ctx); printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(cinv, ctx);
        fq_poly_clear(d, ctx);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
