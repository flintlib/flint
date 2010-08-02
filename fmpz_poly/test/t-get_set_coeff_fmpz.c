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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    printf("get/set_coeff_fmpz....");
    fflush(stdout);

    fmpz_randinit();

    for (i = 0; i < 1000; i++)
    {
        fmpz_poly_t a;
        fmpz_t x1, x2;
        long coeff, len;

        fmpz_poly_init(a);
        fmpz_init(x1);
        fmpz_init(x2);
        len = n_randint(100) + 1;

        for (j = 0; j < 1000; j++)
        {
            fmpz_randtest(x1, 200);
            coeff = n_randint(len);
            fmpz_poly_set_coeff_fmpz(a, coeff, x1);
            fmpz_poly_get_coeff_fmpz(x2, a, coeff);

            result = (fmpz_equal(x1, x2));
            if (!result)
            {
                printf("FAIL: x1 = ");
                fmpz_print(x1);
                printf(", x2 = ");
                fmpz_print(x2);
                printf(", coeff = %ld, length = %ld\n", coeff, len);
                abort();
            }
        }

        fmpz_clear(x1);
        fmpz_clear(x2);
        fmpz_poly_clear(a);
    }

    fmpz_randclear();
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
