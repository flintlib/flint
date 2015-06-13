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

    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"

int main()
{

    fmpz_t a, b, mmin;
    int i, j, k;
    mp_ptr mpna;
    FLINT_TEST_INIT(state);

    flint_printf("get_mpn....");
    fflush(stdout);

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(mmin);


    for (i = 0; i < 100; i++)
    {
        fmpz_set_ui(mmin, 1);
        fmpz_mul_2exp(mmin, mmin, i * 64);

        for (j = 0; j < 100; j++)
        {
            fmpz_set_ui(b, 0);

            k = n_randint(state, 64);
            k += 1;                     /* 1 <= k <= 64 */
            k += (i * 64);              /* 2^(i*64) + 1 <= k <= 2^((i + 1)*64) */

            fmpz_randtest_unsigned(a, state, k);
            fmpz_add(a, a, mmin);

            k = fmpz_get_mpn(&mpna, a);

            while (k)
            {
                fmpz_add_ui(b, b, mpna[k - 1]);

                if (k - 1)
                    fmpz_mul_2exp(b, b, 64);

                k -= 1;
            }

            if (fmpz_cmp(a, b))
            {
                printf("conversion failed.\nn : ");
                fmpz_print(a);
                printf("\nconverted value : ");
                fmpz_print(b);
                printf("\n");
                abort();
            }

            flint_free(mpna);
        }
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(mmin);
    
    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
