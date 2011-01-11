
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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    flint_rand_t state;
    mpq_t n1, n2;

    printf("get/set_coeff_mpq....");
    fflush(stdout);

    mpq_init(n1);
    mpq_init(n2);

    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        fmpq_poly_t a;
        fmpz_t xnum, xden;
        long coeff, len;

        fmpq_poly_init(a);
        fmpz_init(xnum);
        fmpz_init(xden);
        len = (long) (n_randint(100, state) + 1);

        for (j = 0; j < 100; j++)
        {
            fmpz_randtest(xnum, state, 200);
            fmpz_randtest_not_zero(xden, state, 200);
            fmpz_get_mpz(mpq_numref(n1), xnum);
            fmpz_get_mpz(mpq_denref(n1), xden);
            mpq_canonicalize(n1);
            coeff = (long) n_randint(len, state);
            fmpq_poly_set_coeff_mpq(a, coeff, n1);
            fmpq_poly_get_coeff_mpq(n2, a, coeff);

            result = (mpq_equal(n1, n2));
            if (!result)
            {
                gmp_printf
                    ("FAIL: n1 = %Qd, n2 = %Qd, coeff = %ld, length = %ld\n",
                     n1, n2, coeff, len);
                abort();
            }
        }

        fmpz_clear(xnum);
        fmpz_clear(xden);
        fmpq_poly_clear(a);
    }

    mpq_clear(n1);
    mpq_clear(n2);

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
