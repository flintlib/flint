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
    fmpz_t prime1, prime2, prime3, prime4, primeprod, fac, modval, maxa, maxy, a, y;
    int i, j, k, fails;
    FLINT_TEST_INIT(state);

    fmpz_init(prime1);
    fmpz_init(prime2);
    fmpz_init(prime3);
    fmpz_init(prime4);
    fmpz_init(primeprod);
    fmpz_init(fac);
    fmpz_init(modval);
    fmpz_init(a);
    fmpz_init(y);
    fmpz_init(maxa);
    fmpz_init(maxy);

    fails = 0;

    flint_printf("pollard_brent_single....");
    fflush(stdout);

    for (i = 5; i < 36; i += 5)
    {
        for (j = 0; j < 10 * flint_test_multiplier(); j++)
        {

            fmpz_set_ui(prime1, n_randprime(state, i, 1));
            fmpz_set_ui(prime2, n_randprime(state, i, 1));
            fmpz_set_ui(prime3, n_randprime(state, i, 1));
            fmpz_set_ui(prime4, n_randprime(state, i, 1));

            fmpz_mul(prime1, prime1, prime2);
            fmpz_mul(prime3, prime3, prime4);
            fmpz_mul(primeprod, prime1, prime3);

            /* Assigning random values of y and a */
            fmpz_sub_ui(maxa, primeprod, 3);    
            fmpz_randm(a, state, maxa);  
            fmpz_add_ui(a, a, 1);               /* 1 <= a <= n - 3 */
            fmpz_sub_ui(maxa, primeprod, 1);               
            fmpz_randm(y, state, primeprod);
            fmpz_add_ui(y, y, 1);               /* 1 <= y <= n - 1 */

            k = fmpz_factor_pollard_brent_single(fac, primeprod, y, a, UWORD(1) << i);

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
        printf("FAIL : Pollard Rho failed too many times (%d times)\n", fails);
        abort();
    }

    FLINT_TEST_CLEANUP(state);
    fmpz_clear(prime1);
    fmpz_clear(prime2);
    fmpz_clear(prime3);
    fmpz_clear(prime4);
    fmpz_clear(primeprod);
    fmpz_clear(fac);
    fmpz_clear(modval);

    flint_printf("PASS\n");
    return 0;
}
