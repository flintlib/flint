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

#undef ulong /* prevent clash with standard library */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t flint_pseudosquares[] = {17, 73, 241, 1009, 2641, 8089, 18001, 
          53881, 87481, 117049, 515761, 1083289, 3206641, 3818929, 9257329, 
          22000801, 48473881, 48473881, 175244281, 427733329, 427733329, 
          898716289u, 2805544681u, 2805544681u, 2805544681u
#ifndef FLINT64          
          };
#else
          , 10310263441u, 23616331489u, 85157610409u, 85157610409u, 
          196265095009u, 196265095009u, 2871842842801u, 2871842842801u, 
          2871842842801u, 26250887023729u, 26250887023729u, 112434732901969u, 
          112434732901969u, 112434732901969u, 178936222537081u, 
          178936222537081u, 696161110209049u, 696161110209049u,
          2854909648103881u, 6450045516630769u, 6450045516630769u, 
          11641399247947921u, 11641399247947921u, 190621428905186449u,
          196640248121928601u, 712624335095093521u, 1773855791877850321u };
#endif

#if FLINT64
#define FLINT_NUM_PSEUDOSQUARES 52
#else
#define FLINT_NUM_PSEUDOSQUARES 25
#endif

int n_is_prime_pseudosquare(mp_limb_t n)
{
    unsigned int i, j, m1;
    mp_limb_t p, B, NB, exp, mod8;

    if (n < 2UL) return 0;

    if ((n & 1UL) == 0UL)
    {
        return (n == 2UL);
    }

    n_compute_primes(FLINT_PSEUDOSQUARES_CUTOFF);

    for (i = 0; i < FLINT_PSEUDOSQUARES_CUTOFF; i++)
    {
        double ppre;
        p = flint_primes[i];
        if (p*p > n) return 1;
        ppre = flint_prime_inverses[i];
        if (!n_mod2_precomp(n, p, ppre)) return 0;
    }

    B  = flint_primes[FLINT_PSEUDOSQUARES_CUTOFF];
    NB = (n - 1)/B + 1;
    m1 = 0;

    for (i = 0; i < FLINT_NUM_PSEUDOSQUARES; i++)
    {
        if (flint_pseudosquares[i] > NB) break;
    }

    exp = (n - 1)/2;

    for (j = 0; j <= i; j++)
    {
        mp_limb_t mod = n_powmod2(flint_primes[j], exp, n);
        if ((mod != 1UL) && (mod != n - 1)) return 0;
        if (mod == n - 1) m1 = 1;
    }

    mod8 = n % 8;

    if ((mod8 == 3) || (mod8 == 7)) return 1;

    if (mod8 == 5)
    {
        mp_limb_t mod = n_powmod2(2UL, exp, n);
        if (mod == n - 1) return 1;
        printf("Whoah, %lu is a probable prime, but not prime, please report!!\n", n);
        abort();
    }
    else
    {
        if (m1) return 1;
        for (j = i + 1; j < FLINT_NUM_PSEUDOSQUARES + 1; j++)
        {
            mp_limb_t mod = n_powmod2(flint_primes[j], exp, n);
            if (mod == n - 1) return 1;
            if (mod != 1)
            {
                printf("Whoah, %lu is a probable prime, but not prime, please report!!\n", n);
                abort();
            }
        }
        printf("Whoah, %lu is a probable prime, but not prime, please report!!\n", n);
        abort();
    }
}
