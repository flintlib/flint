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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int main(void)
{
    fmpz_t prime1, prime2, primeprod, fac, modval;
    int i, j, k, fails;

    fmpz_init(prime1);
    fmpz_init(prime2);
    fmpz_init(primeprod);
    fmpz_init(fac);
    fmpz_init(modval);
    FLINT_TEST_INIT(state);

    fails = 0;

    flint_printf("ecm....");
    fflush(stdout);

    for (i = 35; i <= 50; i += 5)
    {
        for (j = 0; j < flint_test_multiplier(); j++)
        {
            fmpz_set_ui(prime1, n_randprime(state, i, 1));
            fmpz_set_ui(prime2, n_randprime(state, i, 1));

            fmpz_mul(primeprod, prime1, prime2);

            k = fmpz_factor_ecm(fac, i << 2, 2000, 50000, state, primeprod);

            if (k == 0)
                fails += 1;
            else
            {
                fmpz_mod(modval, primeprod, fac);
                k = fmpz_cmp_ui(modval, 0);
                if (k != 0)
                {
                    printf("FAIL : Wrong factor calculated\n");
                    printf("n : ");
                    fmpz_print(primeprod);
                    printf(" factor calculated : ");
                    fmpz_print(fac);
                    abort();
                }
            }
        }
    }

    if (fails > flint_test_multiplier())
    {
        printf("FAIL : ECM failed too many times (%d times)\n", fails);
        abort();
    }

    FLINT_TEST_CLEANUP(state);
    fmpz_clear(prime1);
    fmpz_clear(prime2);
    fmpz_clear(primeprod);
    fmpz_clear(fac);
    fmpz_clear(modval);

    flint_printf("PASS\n");
    return 0;
}
