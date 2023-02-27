/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2014 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>
#include <stdio.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    fmpz_mod_ctx_t ctx;
    FLINT_TEST_INIT(state);
    
    flint_printf("compose_mod_brent_kung_vec_preinv....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, ainv, b, c;
        fmpz_t p;
        slong l, j, k;
        fmpz_mod_poly_struct * pow, * res;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(ainv, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);

        fmpz_mod_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1, ctx);
        l = 1 + n_randint(state, 20);
        k = 1 + n_randint(state, l);

        fmpz_mod_poly_reverse(ainv, a, a->length, ctx);
        fmpz_mod_poly_inv_series_newton(ainv, ainv, a->length, ctx);
        pow = (fmpz_mod_poly_struct *) flint_malloc((l + k)*sizeof(fmpz_mod_poly_struct));
        res = pow + l;

        fmpz_mod_poly_rem(b, b, a, ctx);
        for (j = 0; j < l; j++)
        {
            fmpz_mod_poly_init(pow + j, ctx);
            fmpz_mod_poly_randtest(pow + j, state, n_randint(state, 20) + 1, ctx);
            fmpz_mod_poly_rem(pow + j, pow + j, a, ctx);
        }

        for (j = 0; j < k; j++)
        {
            fmpz_mod_poly_init(res + j, ctx);
        }

        fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(res, pow, l, k, b, a, ainv, ctx);

        for (j = 0; j < k; j++)
        {
            fmpz_mod_poly_compose_mod(c, pow + j, b, a, ctx);
            if (!fmpz_mod_poly_equal(res + j, c, ctx))
            {
                flint_printf("FAIL (composition):\n");
                flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx); flint_printf("\n");
                flint_printf("res:\n"); fmpz_mod_poly_print(res + j, ctx); flint_printf("\n");
                flint_printf("pow:\n"); fmpz_mod_poly_print(pow + j, ctx); flint_printf("\n");
                flint_printf("b:\n"); fmpz_mod_poly_print(b, ctx); flint_printf("\n");
                flint_printf("c:\n"); fmpz_mod_poly_print(c, ctx); flint_printf("\n");
                flint_printf("j: %wd\n", j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(ainv, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        for (j = 0; j < l; j++)
            fmpz_mod_poly_clear(pow + j, ctx);
        for (j = 0; j < k; j++)
            fmpz_mod_poly_clear(res + j, ctx);
        flint_free(pow);
    }

    fmpz_mod_ctx_clear(ctx);
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
