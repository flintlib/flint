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

    Copyright (C) 2011 Sebastian Pancratz

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

    printf("inv_series... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 1000; i++)
    {
        padic_poly_t a, b, c;
        long n;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_VAL_UNIT);

        padic_poly_init(a);
        padic_poly_init(b);
        padic_poly_init(c);

        padic_poly_randtest(a, state, n_randint(state, 100) + 1, ctx);
        fmpz_remove(a->coeffs, a->coeffs, p);
        if (fmpz_is_zero(a->coeffs))
        {
            fmpz_randtest_not_zero(a->coeffs, state, 20);
            fmpz_remove(a->coeffs, a->coeffs, p);
            padic_poly_reduce(a, ctx);
        }

        padic_poly_set(b, a);
        n = n_randint(state, 100) + 1;

        padic_poly_inv_series(c, b, n, ctx);
        padic_poly_inv_series(b, b, n, ctx);

        result = (padic_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), padic_poly_print(a, ctx), printf("\n\n");
            printf("b = "), padic_poly_print(b, ctx), printf("\n\n");
            printf("c = "), padic_poly_print(c, ctx), printf("\n\n");
            abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /*
        Check correctness: 

        If ord_p(a) = v then we can compute b = a^{-1} mod p^N 
        and we will have a b = 1 mod p^{N-|v|}.  Thus, require 
        that N - |v| > 0.
     */
    for (i = 0; i < 1000; i++)
    {
        padic_ctx_t ctx2;
        padic_poly_t a, b, c;
        long n;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_VAL_UNIT);

        padic_poly_init(a);
        padic_poly_init(b);
        padic_poly_init(c);

        {
            long i, len = n_randint(state, 10) + 1;
            int alloc;
            fmpz_t pow;

            padic_poly_fit_length(a, len);
            _padic_poly_set_length(a, len);
            a->val = n_randint(state, N);
            if (n_randint(state, 2))
                a->val = - a->val;

            alloc = _padic_ctx_pow_ui(pow, N - a->val, ctx);

            for (i = 0; i < len; i++)
                fmpz_randm(a->coeffs + i, state, pow);
            while (fmpz_is_zero(a->coeffs))
                fmpz_randm(a->coeffs, state, pow);
            fmpz_remove(a->coeffs, a->coeffs, p);
            _padic_poly_normalise(a);

            if (alloc)
                fmpz_clear(pow);
        }

        n = n_randint(state, 100) + 1;

        padic_ctx_init(ctx2, p, N - FLINT_ABS(a->val), PADIC_VAL_UNIT);

        padic_poly_inv_series(b, a, n, ctx);
        padic_poly_mul(c, a, b, ctx2);
        padic_poly_truncate(c, n, p);

        result = (padic_poly_is_one(c, ctx2));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), padic_poly_print(a, ctx), printf("\n\n");
            printf("b = "), padic_poly_print(b, ctx), printf("\n\n");
            printf("c = "), padic_poly_print(c, ctx2), printf("\n\n");
            printf("p = %ld\n", *p);
            printf("N = %ld\n", N);
            abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
