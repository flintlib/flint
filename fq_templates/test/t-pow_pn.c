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

    flint_printf("pow_pn...");
    fflush(stdout);

    /* No aliasing */
    for (i = 0; i < 50; i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, ground_mat_t) ap_mat;
        TEMPLATE(T, t) a, res1, res2;
        ulong n;
        fmpz_t exp;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(res1, ctx);
        TEMPLATE(T, init)(res2, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);

        fmpz_init_set(exp, TEMPLATE(T, ctx_prime)(ctx));
        n = n_randint(state, 3) + 1;
        fmpz_pow_ui(exp, exp, n);

        TEMPLATE(T, pow)(res1, a, exp, ctx);

        TEMPLATE(T, pow_pn_init_precomp_matrix)(ap_mat, ctx);
        TEMPLATE(T, pow_pn_precompute_matrix_ui)(ap_mat, n, ctx);
        TEMPLATE(T, pow_pn_precomp)(res2, a, ap_mat, ctx);

        result = (TEMPLATE(T, equal)(res1, res2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, ctx_print)(ctx);
            flint_printf("a:\n"); TEMPLATE(T, print_pretty)(a, ctx);
            flint_printf("\n\n");
            flint_printf("res1:\n"); TEMPLATE(T, print_pretty)(res1, ctx);
            flint_printf("\n\n");
            flint_printf("res2:\n"); TEMPLATE(T, print_pretty)(res2, ctx);
            flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(res1, ctx);
        TEMPLATE(T, clear)(res2, ctx);
        TEMPLATE(T, pow_pn_clear_precomp_matrix)(ap_mat, ctx);
        fmpz_clear(exp);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}


#endif
