/*
    Copyright (C) 2010 Sebastian Pancratz
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
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    ulong cflags = UWORD(0);

    mpq_t n1, n2;

    FLINT_TEST_INIT(state);

    flint_printf("get/set_coeff_mpz....");
    fflush(stdout);

    mpq_init(n1);
    mpq_init(n2);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a;
        fmpz_t x1, x2;
        slong coeff, len;

        fmpq_poly_init(a);
        fmpz_init(x1);
        fmpz_init(x2);
        len = (slong) (n_randint(state, 100) + 1);

        for (j = 0; j < 100; j++)
        {
            fmpz_randtest(x1, state, 200);
            fmpz_get_mpz(mpq_numref(n1), x1);
            flint_mpz_set_si(mpq_denref(n1), 1);
            coeff = (slong) n_randint(state, len);
            fmpq_poly_set_coeff_mpz(a, coeff, mpq_numref(n1));
            fmpq_poly_get_coeff_mpq(n2, a, coeff);

            result = (mpq_equal(n1, n2));
            if (!result)
            {
                flint_printf("FAIL:\n\n");
                flint_printf("a     = "), fmpq_poly_debug(a), flint_printf("\n\n");
                flint_printf("coeff = %wd\n\n", coeff);
                flint_printf("len   = %wd\n\n", len);
                flint_printf("cflags = %wu\n\n", cflags);
                gmp_printf("n1 = %Qd\n\n", n1);
                gmp_printf("n2 = %Qd\n\n", n2);
                abort();
            }
        }

        fmpz_clear(x1);
        fmpz_clear(x2);
        fmpq_poly_clear(a);
    }

    mpq_clear(n1);
    mpq_clear(n2);

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
