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
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    flint_rand_t state;
    ulong n1;
    mpq_t n1mpq, n2mpq;

    printf("get/set_coeff_ui....");
    fflush(stdout);
    
    flint_randinit(state);

    mpq_init(n1mpq);
    mpq_init(n2mpq);
    for (i = 0; i < 1000; i++)
    {
        fmpq_poly_t a;
        long coeff, len;

        fmpq_poly_init(a);
        len = (long) (n_randint(100, state) + 1);
        fmpq_poly_randtest(a, state, len, 100);

        for (j = 0; j < 1000; j++)
        {
            n1 = n_randtest(state);
            mpz_set_ui(mpq_numref(n1mpq), n1);
            mpz_set_si(mpq_denref(n1mpq), 1);
            coeff = (long) n_randint(len, state);
            fmpq_poly_set_coeff_ui(a, coeff, n1);
            fmpq_poly_get_coeff_mpq(n2mpq, a, coeff);

            result = (mpq_equal(n1mpq, n2mpq));
            if (!result)
            {
                gmp_printf
                    ("FAIL: n1 = %Qd, n2 = %Qd, coeff = %ld, length = %ld\n",
                     n1mpq, n2mpq, coeff, len);
                abort();
            }
        }

        fmpq_poly_clear(a);
    }

    mpq_clear(n1mpq);
    mpq_clear(n2mpq);
    
    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
