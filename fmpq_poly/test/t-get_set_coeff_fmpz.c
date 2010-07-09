
/*============================================================================

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

===============================================================================*/

/****************************************************************************

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main (void)
{
    int result;

    printf ("get/set_coeff_fmpz....");
    fflush (stdout);

    ulong i, j, coeff;

    ulong n1, n2;

    mpq_t n1mpq, n2mpq;

    mpq_init (n1mpq);
    mpq_init (n2mpq);

    fmpz_randinit ();

    for (i = 0; i < 1000UL; i++)
    {
        fmpq_poly_t a;

        fmpz_t x1, x2;

        ulong length;

        fmpq_poly_init (a);
        fmpz_init (x1);
        fmpz_init (x2);
        length = n_randint (100) + 1;

        for (j = 0; j < 1000; j++)
        {
            fmpz_randtest (x1, 200);
            fmpz_get_mpz (mpq_numref (n1mpq), x1);
            mpz_set_si (mpq_denref (n1mpq), 1);
            coeff = n_randint (length);
            fmpq_poly_set_coeff_fmpz (a, coeff, x1);
            fmpq_poly_get_coeff_mpq (n2mpq, a, coeff);

            result = (mpq_equal (n1mpq, n2mpq));
            if (!result)
            {
                printf ("Error: n1 = %s, n2 = %s, coeff = %lu, length = %lu\n",
                        mpq_get_str (NULL, 10, n1mpq), mpq_get_str (NULL, 10,
                                                                    n2mpq),
                        coeff, length);
                abort ();
            }
        }

        fmpz_clear (x1);
        fmpz_clear (x2);
        fmpq_poly_clear (a);
    }

    mpq_clear (n1mpq);
    mpq_clear (n2mpq);

    fmpz_randclear ();

    _fmpz_cleanup ();
    printf ("PASS\n");
    return 0;
}
