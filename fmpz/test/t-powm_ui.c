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

    Copyright (C) 2009 William Hart
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

    printf("powm_ui....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with MPIR */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, m;
        ulong x;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(m);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(c, state, 200);
        fmpz_abs(c, c);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(m, c);
        x = n_randtest(state);

        fmpz_powm_ui(b, a, x, c);
        mpz_powm_ui(e, d, x, m);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(e, f) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, x = %lu, m = %Zd\n", d, e, f, x, m);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(m);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        ulong n;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_randtest(b, state, 200);
        fmpz_randtest_not_zero(c, state, 200);
        fmpz_abs(c, c);
        n = n_randtest(state);

        fmpz_powm_ui(a, b, n, c);
        fmpz_powm_ui(b, b, n, c);

        result = (fmpz_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            printf("n = %lu\n", n);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        ulong n;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_randtest(b, state, 200);
        fmpz_randtest_not_zero(c, state, 200);
        fmpz_abs(c, c);
        n = n_randtest(state);

        fmpz_powm_ui(a, b, n, c);
        fmpz_powm_ui(c, b, n, c);

        result = (fmpz_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            printf("n = %lu\n", n);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }

    /* Check aliasing of a and {b, c} */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        ulong n;

        fmpz_init(a);
        fmpz_init(c);

        fmpz_randtest_not_zero(c, state, 200);
        fmpz_abs(c, c);
        n = n_randtest(state);

        fmpz_powm_ui(a, c, n, c);
        fmpz_powm_ui(c, c, n, c);

        result = (fmpz_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            printf("n = %lu\n", n);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
