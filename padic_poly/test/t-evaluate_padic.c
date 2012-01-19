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

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "long_extras.h"
#include "ulong_extras.h"
#include "padic_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    padic_ctx_t ctx;
    fmpz_t p;
    long N;

    printf("evaluate_padic... ");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with the computation over QQ */
    for (i = 0; i < 2000; i++)
    {
        padic_poly_t f;
        fmpq_poly_t fQQ;
        padic_t a, y, z;
        fmpq_t aQQ, yQQ;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_VAL_UNIT);

        padic_poly_init(f);
        fmpq_poly_init(fQQ);
        _padic_init(a);
        _padic_init(y);
        _padic_init(z);
        fmpq_init(aQQ);
        fmpq_init(yQQ);

        padic_poly_randtest(f, state, n_randint(state, 100), ctx);
        padic_randtest(a, state, ctx);

        padic_poly_get_fmpq_poly(fQQ, f, ctx);
        padic_get_fmpq(aQQ, a, ctx);

        padic_poly_evaluate_padic(y, f, a, ctx);
        fmpq_poly_evaluate_fmpq(yQQ, fQQ, aQQ);
        padic_set_fmpq(z, yQQ, ctx);

        if (padic_val(a) >= 0)
        {
            result = (_padic_equal(y, z));
            if (!result)
            {
                printf("FAIL (cmp with QQ):\n");
                printf("f = "), padic_poly_print(f, ctx), printf("\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n\n");
                printf("y = "), padic_print(y, ctx), printf("\n\n");
                printf("z = "), padic_print(z, ctx), printf("\n\n");
                abort();
            }
        }
        else
        {
            padic_ctx_t ctx2;
            padic_t y2, z2;

            padic_ctx_init(ctx2, ctx->p, 
                           ctx->N + (f->length - 1) * padic_val(a), 
                           PADIC_SERIES);
            _padic_init(y2); 
            _padic_init(z2);

            padic_set(y2, y, ctx2);
            padic_set(z2, z, ctx2);

            result = (_padic_equal(y2, z2));
            if (!result)
            {
                printf("FAIL (cmp with QQ):\n");
                printf("f = "), padic_poly_print(f, ctx), printf("\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n\n");
                printf("y = "), padic_print(y, ctx), printf("\n\n");
                printf("z = "), padic_print(z, ctx), printf("\n\n");
                printf("y2 = "), padic_print(y2, ctx2), printf("\n\n");
                printf("z2 = "), padic_print(z2, ctx2), printf("\n\n");
                abort();
            }

            padic_ctx_clear(ctx2);
            _padic_clear(y2);
            _padic_clear(z2);
        }

        padic_poly_clear(f);
        fmpq_poly_clear(fQQ);
        _padic_clear(a);
        _padic_clear(y);
        _padic_clear(z);
        fmpq_clear(aQQ);
        fmpq_clear(yQQ);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
