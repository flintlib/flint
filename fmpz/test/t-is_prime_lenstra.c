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

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

int
main(void)
{
    int i, result, r1;
    slong num = 0;
    FLINT_TEST_INIT(state);

    flint_printf("is_prime_lenstra....");
    fflush(stdout);

    /* test primes always pass n^k - 1 for k = 3, 4, 5, 6, 8, 10, 12 */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        double logd;
        fmpz_t n, F, nk1;
        ulong limit, num_pk1;
        mp_ptr pk1;
        slong k, p;
        n_primes_t iter;
        fmpz * r;

        fmpz_init(n);
        fmpz_init(nk1);
        fmpz_init(F);

        for (k = 3; k <= 12; k += 1 + (k > 5))
        {
           do {
              fmpz_randbits(n, state, n_randint(state, 100) + 2);
              fmpz_abs(n, n);
           } while (!fmpz_is_probabprime(n) || fmpz_cmp_ui(n, 2) == 0);

           logd = log(fmpz_get_d(n));
           limit = (ulong) (logd*logd*logd/10.0) + 20;
   
           pk1 = _nmod_vec_init(k*(ulong) logd);
           num_pk1 = 0;

           fmpz_pow_ui(nk1, n, k);
           fmpz_sub_ui(nk1, nk1, 1);

           n_primes_init(iter);
 
           /* trial factor n^k - 1 up to limit */
           for (p = n_primes_next(iter); p < limit; p = n_primes_next(iter))
           {
              if (fmpz_divisible_si(nk1, p))
                 pk1[num_pk1++] = p;
           }

           r = _fmpz_vec_init(k);

           r1 = fmpz_is_prime_lenstra(F, r, n, pk1, num_pk1, k);

           result = (r1 == 1 && fmpz_divisible(nk1, F));
           if (!result)
           {
               flint_printf("FAIL:\n");
               printf("k = %ld\n", k);
               printf("F = "); fmpz_print(F); printf("\n");
               printf("n = "); fmpz_print(n); printf("\n");
               abort();
           }

           _fmpz_vec_clear(r, k);
           n_primes_clear(iter);
           _nmod_vec_clear(pk1);

        }

        fmpz_clear(n);
        fmpz_clear(nk1);
        fmpz_clear(F);
    }

    /* test composites never pass n^k - 1 for k = 3, 4, 5, 6, 8, 10, 12 */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        double logd;
        fmpz_t n, F, nk1, a;
        ulong limit, num_pk1;
        mp_ptr pk1;
        slong k, p;
        n_primes_t iter;
        fmpz * r;

        fmpz_init(n);
        fmpz_init(a);
        fmpz_init(nk1);
        fmpz_init(F);

        for (k = 3; k <= 12; k += 1 + (k > 5))
        {
           do {
              fmpz_randbits(n, state, n_randint(state, 50) + 2);
           } while (fmpz_cmp_ui(n, 2) < 0 || fmpz_is_even(n));
           do {
              fmpz_randbits(a, state, n_randint(state, 50) + 2);
           } while (fmpz_cmp_ui(a, 2) < 0 || fmpz_is_even(a));
        
           fmpz_mul(n, n, a);

           logd = log(fmpz_get_d(n));
           limit = (ulong) (logd*logd*logd/10.0) + 20;
   
           pk1 = _nmod_vec_init(k*(ulong) logd);
           num_pk1 = 0;

           fmpz_pow_ui(nk1, n, k);
           fmpz_sub_ui(nk1, nk1, 1);

           n_primes_init(iter);
 
           for (p = n_primes_next(iter); p < limit; p = n_primes_next(iter))
           {
              if (fmpz_divisible_si(nk1, p))
                 pk1[num_pk1++] = p;
           }

           r = _fmpz_vec_init(k);

           r1 = fmpz_is_prime_lenstra(F, r, n, pk1, num_pk1, k);

           if (r1 == 1)
              num++;

           result = num < 3; /* almost infinitesimally few pseudoprimes in practice */
           if (!result)
           {
               flint_printf("FAIL:\n");
               printf("Too many pseudoprimes\n");
               printf("k = %ld\n", k);
               printf("n = "); fmpz_print(n); printf("\n");
               printf("F = "); fmpz_print(F); printf("\n");
               abort();
           }

           _fmpz_vec_clear(r, k);
           n_primes_clear(iter);
           _nmod_vec_clear(pk1);
        }

        fmpz_clear(n);
        fmpz_clear(a);
        fmpz_clear(nk1);
        fmpz_clear(F);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
