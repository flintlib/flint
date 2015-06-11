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

/* Test code, factor_ecm_one vs factor_ecm_two */

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "time.h"

int main(int argc, char const *argv[])
{
    fmpz_t prime1, prime2, prod, mod, f;
    int i, j, k, fails1, fails2, res;
    float time1, time2, temp;
    clock_t begin, end;

    fmpz_init(f);
    fmpz_init(prime1);
    fmpz_init(prime2);
    fmpz_init(prod);
    fmpz_init(mod);

    FLINT_TEST_INIT(state);

    printf("\n\n************** TIMING FACTOR_ECM_ONE VS FACTOR_ECM_TWO (100 semiprimes each) ************** \n\n");
    printf("bits1        bits2         time1         fails1         time2         fails2         ratio\n");

    for(i = 30; i <= 61; i += 5)
    {
        for (j = i; j <= 61; j += 5)
        {
            time1 = 0.0;
            time2 = 0.0;

            int fails1 = 0;
            int fails2 = 0;

            for (k = 0; k < 100; k++)
            {   
                fmpz_set_ui(prime1, n_randprime(state, i, 1));
                fmpz_set_ui(prime2, n_randprime(state, j, 1));
                fmpz_mul(prod, prime1, prime2);

                /************************ TIME ecm_factor_one ************************/

                begin = clock();

                res = fmpz_factor_ecm_one(f, (i + j) << 2, 5000, 50000, state, prod);

                end = clock();
                temp = (double)(end - begin) / CLOCKS_PER_SEC;
                time1 += temp;

                fmpz_mod(mod, prod, f);

                if (res)
                {
                    if (fmpz_cmp_ui(mod, 0) || !fmpz_cmp(f, prod) || !fmpz_cmp_ui(f, 1))
                    {
                        printf("WRONG ANSWER\n");
                        abort();
                    }
                }
                else
                    fails1 += 1;

                /************************ TIME ecm_factor_two ************************/

                begin = clock();

                if (!fmpz_factor_ecm_two(f, (i + j) << 2, 5000, 50000, state, prod))
                        fails2 += 1;

                end = clock();
                temp = (double)(end - begin) / CLOCKS_PER_SEC;
                time2 += temp;

                fmpz_mod(mod, prod, f);

                if (res)
                {
                    if (fmpz_cmp_ui(mod, 0) || !fmpz_cmp(f, prod) || !fmpz_cmp_ui(f, 1))
                    {
                        printf("WRONG ANSWER\n");
                        abort();
                    }
                }
                else
                    fails2 += 1;
            }
            printf("  %d          %d          %f          %d          %f          %d          %f\n", i, j, time1, fails1, time2, fails2, time1/time2);

        }
    }

    return 0;
}
