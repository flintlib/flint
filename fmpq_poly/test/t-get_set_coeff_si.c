/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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
#include "fmpq_poly.h"
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    ulong cflags = UWORD(0);

    slong n;
    mpq_t n_mpq;
    FLINT_TEST_INIT(state);
    

    flint_printf("get/set_coeff_si....");
    fflush(stdout);

    mpq_init(n_mpq);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a;
        slong coeff, len;

        fmpq_poly_init(a);
        len = (slong) n_randint(state, 100) + 1;

        for (j = 0; j < 1000; j++)
        {
            n = z_randtest(state);
            coeff = n_randint(state, len);
            fmpq_poly_set_coeff_si(a, coeff, n);
            fmpq_poly_get_coeff_mpq(n_mpq, a, coeff);

            cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
            result = (flint_mpz_cmp_ui(mpq_denref(n_mpq), 1) == 0 
                   && flint_mpz_cmp_si(mpq_numref(n_mpq), n) == 0
                   && !cflags);
            if (!result)
            {
                flint_printf("FAIL:\n");
                flint_printf("a      = "), fmpq_poly_debug(a), flint_printf("\n");
                flint_printf("len    = %wd\n", len);
                flint_printf("coeff  = %wd\n", coeff);
                flint_printf("cflags = %wu\n", cflags);
                flint_printf("n      = %wd\n", n);
                gmp_printf("n_mpq  = %Qd\n", n_mpq);
                fflush(stdout);
                flint_abort();
            }
        }
        fmpq_poly_clear(a);
    }

    
    mpq_clear(n_mpq);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
