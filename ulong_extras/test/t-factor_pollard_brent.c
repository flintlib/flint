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
#include "ulong_extras.h"

int main(void)
{
    mp_limb_t prime1, prime2, primeprod, fac, modval;
    int i, j, k, l, fails;

    FLINT_TEST_INIT(state);

    fails = 0;

    flint_printf("pollard_brent....");
    fflush(stdout);

    for (l = 5; l < 26; l += 5)    
    {
        for (i = l; i < 26; i += 5)
        {
            for (j = 0; j < 10 * flint_test_multiplier(); j++)
            {
                prime1 = n_randtest_bits(state, l);
                prime2 = n_randtest_bits(state, i);
                primeprod = prime1 * prime2;

                k = n_factor_pollard_brent(&fac, state, primeprod, 5, 2500);

                if (k == 0)
                    fails += 1;
                else
                {
                    modval = primeprod % fac;
                    if (modval != 0)
                    {
                        flint_printf("FAIL : Wrong factor calculated\n");
                        flint_printf("n : %wu\n", primeprod);
                        flint_printf("Factor calculated: %wn\n", fac);
                        abort();
                    }
                }
            }
        }
    }

    if (fails > flint_test_multiplier())
    {
        printf("FAIL : Pollard Rho - Brent failed too many times (%d times)\n", fails);
        abort();
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
