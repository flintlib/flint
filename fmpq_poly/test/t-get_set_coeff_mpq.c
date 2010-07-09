
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
main (void)
{
    int result;

    printf ("get/set_coeff_mpq....");
    fflush (stdout);

    ulong i, j, coeff;

    mpq_t n1mpq, n2mpq;

    mpq_init (n1mpq);
    mpq_init (n2mpq);

    fmpz_randinit ();

    for (i = 0; i < 100UL; i++)
    {
        fmpq_poly_t a;

        fmpz_t xnum, xden;

        ulong length;

        fmpq_poly_init (a);
        fmpz_init (xnum);
        fmpz_init (xden);
        length = n_randint (100) + 1;

        for (j = 0; j < 100UL; j++)
        {
            fmpz_randtest (xnum, 200);
            fmpz_randtest_not_zero (xden, 200);
            fmpz_get_mpz (mpq_numref (n1mpq), xnum);
            fmpz_get_mpz (mpq_denref (n1mpq), xden);
            mpq_canonicalize (n1mpq);
            coeff = n_randint (length);
            fmpq_poly_set_coeff_mpq (a, coeff, n1mpq);
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

        fmpz_clear (xnum);
        fmpz_clear (xden);
        fmpq_poly_clear (a);
    }

    fmpz_randclear ();

    mpq_clear (n1mpq);
    mpq_clear (n2mpq);

    _fmpz_cleanup ();
    printf ("PASS\n");
    return 0;
}
