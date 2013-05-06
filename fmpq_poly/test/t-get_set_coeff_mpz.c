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
    Copyright (C) 2009 William Hart

******************************************************************************/

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
    flint_rand_t state;
    ulong cflags = 0UL;

    mpq_t n1, n2;

    printf("get/set_coeff_mpz....");
    fflush(stdout);

    mpq_init(n1);
    mpq_init(n2);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a;
        fmpz_t x1, x2;
        long coeff, len;

        fmpq_poly_init(a);
        fmpz_init(x1);
        fmpz_init(x2);
        len = (long) (n_randint(state, 100) + 1);

        for (j = 0; j < 100; j++)
        {
            fmpz_randtest(x1, state, 200);
            fmpz_get_mpz(mpq_numref(n1), x1);
            mpz_set_si(mpq_denref(n1), 1);
            coeff = (long) n_randint(state, len);
            fmpq_poly_set_coeff_mpz(a, coeff, mpq_numref(n1));
            fmpq_poly_get_coeff_mpq(n2, a, coeff);

            result = (mpq_equal(n1, n2));
            if (!result)
            {
                printf("FAIL:\n\n");
                printf("a     = "), fmpq_poly_debug(a), printf("\n\n");
                printf("coeff = %ld\n\n", coeff);
                printf("len   = %ld\n\n", len);
                printf("cflags = %lu\n\n", cflags);
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

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
