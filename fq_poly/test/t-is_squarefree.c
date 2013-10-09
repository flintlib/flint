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
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "fq_poly.h"

int
main(void)
{
    int iter;
    flint_rand_t state;
    flint_randinit(state);

    flint_printf("is_squarefree....");
    fflush(stdout);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        fq_ctx_t ctx;
        fq_poly_t poly, Q, R, t;
        fmpz_t x;
        slong i, num_factors, exp, max_exp;
        int v, result;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(poly);
        fq_poly_init(t);
        fq_poly_init(Q);
        fq_poly_init(R);

        fmpz_init(x);
        fmpz_randtest_mod(x, state, fq_ctx_prime(ctx));

        fq_poly_set_coeff_fmpz(poly, 0, x, ctx);
        num_factors = n_randint(state, 5);

        max_exp = 0;
        for (i = 0; i < num_factors; i++)
        {
            do {
                fq_poly_randtest(t, state, n_randint(state, 10), ctx);
            } while (!fq_poly_is_irreducible(t, ctx) ||
                    (fq_poly_length(t) < 2));

            exp = n_randint(state, 4) + 1;
            if (n_randint(state, 2) == 0)
                exp = 1;

            fq_poly_divrem(Q, R, poly, t, ctx);
            if (!fq_poly_is_zero(R))
            {
                fq_poly_pow(t, t, exp, ctx);
                fq_poly_mul(poly, poly, t, ctx);
                max_exp = FLINT_MAX(exp, max_exp);
            }
        }

        v = fq_poly_is_squarefree(poly, ctx);

        if (v == 1)
            result = (max_exp <= 1 && !fq_poly_is_zero(poly));
        else
            result = (max_exp > 1 || fq_poly_is_zero(poly));

        if (!result)
        {
            flint_printf("FAIL: ");
            fq_ctx_print(ctx);
            flint_printf(" %ld, %d\n", max_exp, v);
            fq_poly_print(poly, ctx); flint_printf("\n");
            abort();
        }

        fq_poly_clear(poly);
        fq_poly_clear(t);
        fq_poly_clear(Q);
        fq_poly_clear(R);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return 0;
}
