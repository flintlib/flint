/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
    FLINT_TEST_INIT(state);

    flint_printf("scalar_fdiv_mpz....");
    fflush(stdout);

    

    /* Compare with fmpz_poly_scalar_fdiv_fmpz */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        fmpz_t n;
        mpz_t n1;

        fmpz_init(n);
        mpz_init(n1);
        do {
            fmpz_randtest(n, state, 200);
        } while (fmpz_is_zero(n));
        fmpz_get_mpz(n1, n);

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_scalar_mul_fmpz(a, a, n);

        fmpz_poly_scalar_fdiv_fmpz(b, a, n);
        fmpz_poly_scalar_fdiv_mpz(c, a, n1);

        result = (fmpz_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            abort();
        }

        /* aliasing */
        fmpz_poly_scalar_fdiv_mpz(a, a, n1);
        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            abort();
        }

        mpz_clear(n1);
        fmpz_clear(n);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
