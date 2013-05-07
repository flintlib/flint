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

    printf("remove....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with MPIR, random input */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, g;
        long x, y;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest_not_zero(a, state, 200);
        do {
            fmpz_randtest_not_zero(b, state, 200);
            fmpz_abs(b, b);
        } while (fmpz_is_one(b));
        

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        x = fmpz_remove(c, a, b);
        y = mpz_remove(f, d, e);

        fmpz_get_mpz(g, c);

        result = ((x == y) && (mpz_cmp(f, g) == 0));

        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd\n", d, e, f, g);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Compare with MPIR, random input but ensure that factors exist */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, pow;
        mpz_t d, e, f, g;
        long x, y;
        ulong n;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(pow);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest_not_zero(a, state, 200);
        do {
            fmpz_randtest_not_zero(b, state, 200);
            fmpz_abs(b, b);
        } while (fmpz_is_one(b));

        n = n_randint(state, 10);
        fmpz_pow_ui(pow, b, n);
        fmpz_mul(a, a, pow);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        x = fmpz_remove(c, a, b);
        y = mpz_remove(f, d, e);

        fmpz_get_mpz(g, c);

        result = ((x == y) && (x >= n) && (mpz_cmp(f, g) == 0));

        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd\n", d, e, f, g);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(pow);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        long x;

        fmpz_init(a);
        fmpz_init(c);

        do {
            fmpz_randtest_not_zero(a, state, 200);
            fmpz_abs(a, a);
        } while (fmpz_is_one(a));

        x = fmpz_remove(c, a, a);

        result = ((x == 1) && (fmpz_cmp_ui(c, 1) == 0));

        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(a), printf("\n");
            fmpz_print(c), printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        long x, y;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_randtest_not_zero(a, state, 200);
        do {
            fmpz_randtest_not_zero(b, state, 200);
            fmpz_abs(b, b);
        } while (fmpz_is_one(b));

        x = fmpz_remove(c, a, b);
        y = fmpz_remove(a, a, b);

        result = ((x == y) && fmpz_equal(a, c));

        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(a), printf("\n");
            fmpz_print(b), printf("\n");
            fmpz_print(c), printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        long x, y;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_randtest_not_zero(a, state, 200);
        do {
            fmpz_randtest_not_zero(b, state, 200);
            fmpz_abs(b, b);
        } while (fmpz_is_one(b));

        x = fmpz_remove(c, a, b);
        y = fmpz_remove(b, a, b);

        result = ((x == y) && fmpz_equal(b, c));

        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(a), printf("\n");
            fmpz_print(b), printf("\n");
            fmpz_print(c), printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
