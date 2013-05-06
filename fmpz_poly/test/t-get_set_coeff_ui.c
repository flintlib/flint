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

Copyright (C) 2009 William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    flint_rand_t state;

    printf("get/set_coeff_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a;
        ulong n1, n2;
        long coeff, len;

        fmpz_poly_init(a);
        len = n_randint(state, 100) + 1;

        for (j = 0; j < 1000; j++)
        {
            n1 = n_randtest(state);
            coeff = n_randint(state, len);
            fmpz_poly_set_coeff_ui(a, coeff, n1);
            n2 = fmpz_poly_get_coeff_ui(a, coeff);

            result = (n1 == n2);
            if (!result)
            {
                printf("FAIL: n1 = %lu, n2 = %lu, coeff = %ld, length = %ld\n",
                       n1, n2, coeff, len);
                abort();
            }
        }

        fmpz_poly_clear(a);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
