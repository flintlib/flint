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
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("clog....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        long l;

        fmpz_init(a);

        while (fmpz_cmp_ui(a, 2) < 0)
            fmpz_randtest(a, state, 200);

        l = fmpz_clog(a, a);

        result = (l == 1);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("l = %ld\n", l);
            abort();
        }

        fmpz_clear(a);
    }

    /* Check correctness */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, x, y;
        long k;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(x);
        fmpz_init(y);

        while (fmpz_cmp_ui(a, 1) < 0)
            fmpz_randtest(a, state, 200);
        while (fmpz_cmp_ui(b, 2) < 0)
            fmpz_randtest(b, state, 200);

        k = fmpz_clog(a, b);  /* p^{k-1} < a <= p^k*/

        if (k > 0)
            fmpz_pow_ui(x, b, k - 1);
        else
            fmpz_zero(x);
        fmpz_pow_ui(y, b, k);

        result = (fmpz_cmp(x, a) < 0 && fmpz_cmp(a, y) <= 0);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("x = "), fmpz_print(x), printf("\n");
            printf("y = "), fmpz_print(y), printf("\n");
            printf("k = %ld\n", k);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(x);
        fmpz_clear(y);
    }

    /* Check correctness:  exact powers */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        long k, l;

        fmpz_init(a);
        fmpz_init(b);

        while (fmpz_cmp_ui(b, 2) < 0)
            fmpz_randtest(b, state, 200);
        l = n_randint(state, 20);
        fmpz_pow_ui(a, b, l);

        k = fmpz_clog(a, b);

        result = (k == l);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("k = %ld\n", k);
            printf("l = %ld\n", l);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

