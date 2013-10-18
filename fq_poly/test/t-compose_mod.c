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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fq_poly.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);
    flint_printf("compose_mod....");
    fflush(stdout);

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, c, d, e;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(d, ctx);
        fq_poly_init(e, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fq_poly_compose_mod(d, a, b, c, ctx);
        fq_poly_compose(e, a, b, ctx);
        fq_poly_rem(e, e, c, ctx);

        if (!fq_poly_equal(d, e, ctx))
        {
            flint_printf("FAIL (composition):\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx); flint_printf("\n");
            flint_printf("b:\n"); fq_poly_print(b, ctx); flint_printf("\n");
            flint_printf("c:\n"); fq_poly_print(c, ctx); flint_printf("\n");
            flint_printf("d:\n"); fq_poly_print(d, ctx); flint_printf("\n");
            flint_printf("e:\n"); fq_poly_print(e, ctx); flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(d, ctx);
        fq_poly_clear(e, ctx);

        fq_ctx_clear(ctx);
    }

    /* Test aliasing of res and a */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, c, d;

        fq_ctx_randtest(ctx, state);
        
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(d, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fq_poly_compose_mod(d, a, b, c, ctx);
        fq_poly_compose_mod(a, a, b, c, ctx);

        if (!fq_poly_equal(d, a, ctx))
        {
            flint_printf("FAIL (aliasing a):\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx); flint_printf("\n");
            flint_printf("b:\n"); fq_poly_print(b, ctx); flint_printf("\n");
            flint_printf("c:\n"); fq_poly_print(c, ctx); flint_printf("\n");
            flint_printf("d:\n"); fq_poly_print(d, ctx); flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(d, ctx);

        fq_ctx_clear(ctx);
    }

    /* Test aliasing of res and b */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, c, d;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(d, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fq_poly_compose_mod(d, a, b, c, ctx);
        fq_poly_compose_mod(b, a, b, c, ctx);

        if (!fq_poly_equal(d, b, ctx))
        {
            flint_printf("FAIL (aliasing b)\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx); flint_printf("\n");
            flint_printf("b:\n"); fq_poly_print(b, ctx); flint_printf("\n");
            flint_printf("c:\n"); fq_poly_print(c, ctx); flint_printf("\n");
            flint_printf("d:\n"); fq_poly_print(d, ctx); flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(d, ctx);

        fq_ctx_clear(ctx);
    }

    /* Test aliasing of res and c */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, c, d;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(d, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fq_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fq_poly_compose_mod(d, a, b, c, ctx);
        fq_poly_compose_mod(c, a, b, c, ctx);

        if (!fq_poly_equal(d, c, ctx))
        {
            flint_printf("FAIL (aliasing c)\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx); flint_printf("\n");
            flint_printf("b:\n"); fq_poly_print(b, ctx); flint_printf("\n");
            flint_printf("c:\n"); fq_poly_print(c, ctx); flint_printf("\n");
            flint_printf("d:\n"); fq_poly_print(d, ctx); flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(d, ctx);
        
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return 0;
}
