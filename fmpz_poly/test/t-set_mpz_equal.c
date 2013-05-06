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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    gmp_randstate_t state;

    printf("set_mpz_equal....");
    fflush(stdout);

    gmp_randinit_default(state);
    gmp_randseed_ui(state, 23);

    /* equal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        mpz_t n;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        mpz_init(n);

        mpz_rrandomb(n, state, 200);
        fmpz_poly_set_mpz(a, n);
        fmpz_poly_set(b, a);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("n = %Zd\n\n", n);
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        mpz_clear(n);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        mpz_t m, n;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        mpz_init(m);
        mpz_init(n);

        mpz_rrandomb(m, state, 200);
        mpz_rrandomb(n, state, 200);
        while (mpz_cmp(m, n) == 0)
            mpz_rrandomb(n, state, 200);
        fmpz_poly_set_mpz(a, m);
        fmpz_poly_set_mpz(b, n);

        result = (!fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("m = %Zd\n\n", m);
            gmp_printf("n = %Zd\n\n", n);
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        mpz_clear(m);
        mpz_clear(n);
    }

    gmp_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
