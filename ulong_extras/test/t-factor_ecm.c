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
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    int i, j, k, result, fails;
    mp_limb_t prime1, prime2, prod, f, mod;
    FLINT_TEST_INIT(state);

    fails = 0;

    flint_printf("factor_ecm....");
    fflush(stdout);

    for (i = 10; i < 64; i += 5)
    {
        for (j = i; j < 64 - i; j += 5)
        {
            for (k = 0; k < flint_test_multiplier(); k++)
            {
                prime1 = n_randprime(state, i, 1);
                prime2 = n_randprime(state, j, 1);
                prod = prime1 * prime2;

                result = n_factor_ecm(&f, (i + j) << 2, 1000, 50000, state, prod);

                if (result)
                {
                    mod = prod % f;
                    if ((mod != 0) || (f == prod) || (f == 1))
                    {
                        flint_printf("WRONG ANSWER from stage %d\n", result);
                        flint_printf("Number : %wu = %wu * %wu\n", prod, prime1, prime2);
                        flint_printf("Factor found : %wu", f);
                        flint_printf("Aborting");
                        abort();
                    }
                }
                else
                    fails += 1;
            }
        }
    }

    if (fails > flint_test_multiplier())
    {
        flint_printf("Too many unsuccessful factorizations, %d\n", fails);
        flint_printf("Aborting\n");
        abort();
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
