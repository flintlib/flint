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
    Copyright (C) 2010 Sebastian Pancratz

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
    long n1;
    mpq_t n1mpq, n2mpq;

    printf("get/set_coeff_si....");
    fflush(stdout);

    mpq_init(n1mpq);
    mpq_init(n2mpq);
    for (i = 0; i < 1000UL; i++)
    {
        fmpq_poly_t a;
        long coeff, len;

        fmpq_poly_init(a);
        len = (long) n_randint(100) + 1;

        for (j = 0; j < 1000; j++)
        {
            n1 = (long) n_randtest();
            coeff = (long) n_randint(len);
            mpq_set_si(n1mpq, n1, 1);
            fmpq_poly_set_coeff_si(a, coeff, n1);
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
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
