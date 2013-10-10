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

    flint_printf("factor_squarefree....");
    fflush(stdout);

    for (iter = 0; iter < 500; iter++)
    {
        int result = 1;
        fq_ctx_t ctx;
        fq_poly_t pol1, poly, quot, rem;
        fq_poly_factor_t res;
        slong exp[5], prod1;
        slong length, i, j, num;

        fq_ctx_randtest(ctx, state);
        
        fq_poly_init(pol1);
        fq_poly_init(poly);
        fq_poly_init(quot);
        fq_poly_init(rem);

        fq_poly_one(pol1, ctx);

        length = n_randint(state, 7) + 2;

        do
        {
            fq_poly_randtest(poly, state, length, ctx);
            fq_poly_make_monic(poly, poly, ctx);
        }
        while ((poly->length != length) || (!fq_poly_is_irreducible(poly, ctx)));
        exp[0] = n_randprime(state, 5, 0);

        prod1 = exp[0];
        for (i = 0; i < exp[0]; i++)
            fq_poly_mul(pol1, pol1, poly, ctx);

        num = n_randint(state, 5) + 1;
        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 7) + 2;
                fq_poly_randtest(poly, state, length, ctx);
                if (poly->length)
                {
                    fq_poly_make_monic(poly, poly, ctx);
                    fq_poly_divrem(quot, rem, pol1, poly, ctx);
                }
            }
            while ((!fq_poly_is_irreducible(poly, ctx)) ||
                   (poly->length != length) || (rem->length == 0));

            do
                exp[i] = n_randprime(state, 5, 0);
            while (prod1 % exp[i] == 0);

            prod1 *= exp[i];
            for (j = 0; j < exp[i]; j++)
                fq_poly_mul(pol1, pol1, poly, ctx);
        }

        fq_poly_factor_init(res, ctx);
        fq_poly_factor_squarefree(res, pol1, ctx);

        result &= (res->num == num);
        if (result)
        {
            ulong prod2 = 1;
            for (i = 0; i < num; i++)
                prod2 *= res->exp[i];
            result &= (prod1 == prod2);
        }

        if (!result)
        {
            flint_printf("Error: exp don't match. Ctx = ");
            fq_ctx_print(ctx);
            flint_printf("\n");
            for (i = 0; i < res->num; i++)
                flint_printf("%ld ", res->exp[i]);
            flint_printf("\n");
            for (i = 0; i < num; i++)
                flint_printf("%ld ", exp[i]);
            flint_printf("\n");
            abort();
        }

        fq_poly_clear(quot);
        fq_poly_clear(rem);
        fq_poly_clear(pol1);
        fq_poly_clear(poly);
        fq_poly_factor_clear(res);
        
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return 0;
}
