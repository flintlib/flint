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

#include "ulong_extras.h"
#include "long_extras.h"
#include "padic.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("randtest... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check randtest() */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t lo, hi, N;
        padic_ctx_t ctx;

        padic_t a;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_randtest(a, state, ctx);

        if (N > 0)
        {
            lo = -((N + 9) / 10);
            hi = N;
        }
        else if (N < 0)
        {
            lo = N - ((-N + 9) / 10);
            hi = N;
        }
        else
        {
            lo = -10;
            hi = 0;
        }

        result = padic_is_zero(a) || (lo <= padic_val(a) && padic_val(a) < hi);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("N = %ld\n", N);
            abort();
        }

        padic_clear(a);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check randtest_not_zero() */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t lo, hi, N;
        padic_ctx_t ctx;

        padic_t a;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_randtest_not_zero(a, state, ctx);

        if (N > 0)
        {
            lo = -((N + 9) / 10);
            hi = N;
        }
        else if (N < 0)
        {
            lo = N - ((-N + 9) / 10);
            hi = N;
        }
        else
        {
            lo = -10;
            hi = 0;
        }

        result = !padic_is_zero(a) && (lo <= padic_val(a) && padic_val(a) < hi);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("N = %ld\n", N);
            abort();
        }

        padic_clear(a);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

