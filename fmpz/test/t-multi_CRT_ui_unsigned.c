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

    Copyright (C) 2008, 2009, William Hart
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"


int main()
{
    int result;
    fmpz_t input, temp;
    mpz_t num1;
    mp_limb_t * output, * output2;
    double primes_per_limb;
    slong i, j, k;
    slong bits;

    slong num_primes;
    mp_limb_t * primes;
    mp_limb_t prime;

    fmpz_comb_t comb;
    fmpz_comb_temp_t comb_temp;

    FLINT_TEST_INIT(state);

    flint_printf("multi_CRT_ui_unsigned....");
    fflush(stdout);

    mpz_init(num1);
    

    result = 1;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        bits = n_randint(state, 300)+1;

        if (FLINT_BITS == 32)
            primes_per_limb = 1.0325;
        else if (FLINT_BITS == 64)
            primes_per_limb = 1.016;

        num_primes = (bits*primes_per_limb)/FLINT_BITS + 1;

        primes = (mp_limb_t *) flint_malloc(num_primes * sizeof(mp_limb_t));
        prime = n_nextprime((UWORD(1) << (FLINT_BITS-1)) - WORD(10000000), 0);

        for (j = 0; j < num_primes; j++)
        {
            primes[j] = prime;
            prime = n_nextprime(prime, 0);
        }

        fmpz_init(input);

        fmpz_randtest(input, state, bits);
        fmpz_abs(input, input);
        fmpz_get_mpz(num1, input);

        output = (mp_limb_t *) flint_malloc(num_primes * sizeof(mp_limb_t));
        output2 = (mp_limb_t *) flint_malloc(num_primes * sizeof(mp_limb_t));

        fmpz_comb_init(comb, primes, num_primes);
        fmpz_comb_temp_init(comb_temp, comb);

        fmpz_multi_mod_ui(output, input, comb, comb_temp);
     
        fmpz_init(temp);

        fmpz_multi_CRT_ui(temp, output, comb, comb_temp, 0);
        result &= fmpz_equal(temp, input);

        fmpz_comb_temp_clear(comb_temp);

        for (k = 0; k < num_primes; k++)
        {
            output2[k] = fmpz_mod_ui(temp, input, primes[k]);
                result &= (output[k] == output2[k]);
        }

        if (!result)
        {
            flint_printf("FAIL: bits = %wd, num_primes = %wd\n", bits, num_primes);
            fmpz_print(temp); flint_printf("\n");
            fmpz_print(input); flint_printf("\n");
            abort();
        }

        fmpz_clear(temp);

        fmpz_comb_clear(comb);
        fmpz_clear(input);

        flint_free(output);
        flint_free(output2);
        flint_free(primes);
    }

    
    mpz_clear(num1);
    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
