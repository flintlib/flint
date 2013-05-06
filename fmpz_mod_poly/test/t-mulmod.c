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

    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("mulmod....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of res and a */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, res, f;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(res, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1);

        fmpz_mod_poly_mulmod(res, a, b, f);
        fmpz_mod_poly_mulmod(a, a, b, f);

        result = (fmpz_mod_poly_equal(res, a));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); fmpz_mod_poly_print(a), printf("\n\n");
            printf("b:\n"); fmpz_mod_poly_print(b), printf("\n\n");
            printf("f:\n"); fmpz_mod_poly_print(f), printf("\n\n");
            printf("res:\n"); fmpz_mod_poly_print(res), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(res);
        fmpz_clear(p);
    }

    /* Check aliasing of res and b */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, f, res;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(res, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1);

        fmpz_mod_poly_mulmod(res, a, b, f);
        fmpz_mod_poly_mulmod(b, a, b, f);

        result = (fmpz_mod_poly_equal(res, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); fmpz_mod_poly_print(a), printf("\n\n");
            printf("b:\n"); fmpz_mod_poly_print(b), printf("\n\n");
            printf("f:\n"); fmpz_mod_poly_print(f), printf("\n\n");
            printf("res:\n"); fmpz_mod_poly_print(res), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(res);
        fmpz_clear(p);
    }

    /* Check aliasing of res and f */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, f, res;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(res, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1);

        fmpz_mod_poly_mulmod(res, a, b, f);
        fmpz_mod_poly_mulmod(f, a, b, f);

        result = (fmpz_mod_poly_equal(res, f));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); fmpz_mod_poly_print(a), printf("\n\n");
            printf("b:\n"); fmpz_mod_poly_print(b), printf("\n\n");
            printf("f:\n"); fmpz_mod_poly_print(f), printf("\n\n");
            printf("res:\n"); fmpz_mod_poly_print(res), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(res);
        fmpz_clear(p);
    }

    /* No aliasing */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, res1, res2, t, f;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1);

        fmpz_mod_poly_init(res1, p);
        fmpz_mod_poly_init(res2, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_mulmod(res1, a, b, f);
        fmpz_mod_poly_mul(res2, a, b);
        fmpz_mod_poly_divrem(t, res2, res2, f);

        result = (fmpz_mod_poly_equal(res1, res2));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); fmpz_mod_poly_print(a), printf("\n\n");
            printf("b:\n"); fmpz_mod_poly_print(b), printf("\n\n");
            printf("f:\n"); fmpz_mod_poly_print(f), printf("\n\n");
            printf("res1:\n"); fmpz_mod_poly_print(res1), printf("\n\n");
            printf("res2:\n"); fmpz_mod_poly_print(res2), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(res1);
        fmpz_mod_poly_clear(res2);
        fmpz_mod_poly_clear(t);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
