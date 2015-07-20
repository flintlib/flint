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

    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#include <stdlib.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

unsigned int nextmod30[] = 
{
   1, 6, 5, 4, 3, 2, 1, 4, 3, 2, 1, 2, 1, 4, 3, 2, 1, 2, 1,
   4, 3, 2, 1, 6, 5, 4, 3, 2, 1, 2
};

unsigned int nextindex[] = 
{
   1, 7, 7, 7, 7, 7, 7, 11, 11, 11, 11, 13, 13, 17, 17, 17, 17, 19, 19,
   23, 23, 23, 23, 29, 29, 29, 29, 29, 29, 1
};

#define NEXTPRIME_PRIMES 54

/* first 64 primes used for modular arithmetic */
static const unsigned short n_modular_primes_tab[64] = {
#if FLINT_BITS == 64
  29, 99, 123, 131, 155, 255, 269, 359, 435, 449, 453, 485, 491, 543, 585,
  599, 753, 849, 879, 885, 903, 995, 1209, 1251, 1311, 1373, 1403, 1485, 1533,
  1535, 1545, 1551, 1575, 1601, 1625, 1655, 1701, 1709, 1845, 1859, 1913,
  1995, 2045, 2219, 2229, 2321, 2363, 2385, 2483, 2499, 2523, 2543, 2613,
  2639, 2679, 2829, 2931, 3089, 3165, 3189, 3245, 3273, 3291, 3341
#else
  11, 45, 65, 95, 129, 135, 165, 209, 219, 221, 239, 245, 281, 303, 345, 351,
  359, 389, 393, 395, 413, 435, 461, 513, 519, 549, 555, 573, 575, 585, 591,
  611, 623, 629, 683, 689, 701, 729, 785, 791, 813, 843, 851, 869, 879, 893,
  905, 921, 953, 963, 965, 969, 993, 1031, 1049, 1073, 1085, 1103, 1143, 1173,
  1203, 1221, 1229, 1271
#endif
};

#define N_MODULUS (UWORD(1) << (FLINT_BITS - 1))

mp_limb_t n_nextprime(mp_limb_t n, int proved)
{
    mp_limb_t * moduli;
    ulong i, index;

    if (n < 7) 
    {
        if (n < 2)
            return 2;
        n++;
        n |= 1;
        return n;  
    }

    if (n >= N_MODULUS && n < N_MODULUS + n_modular_primes_tab[63])
    {
        for (i = 0; i < 64; i++)
            if (N_MODULUS + n_modular_primes_tab[i] > n)
                return N_MODULUS + n_modular_primes_tab[i];
    }

    if (n >= UWORD_MAX_PRIME)
    {
        flint_printf("Exception (n_nextprime). No larger single-limb prime exists.\n");
        abort();
    }

    index = n % 30;
    n += nextmod30[index];
    index = nextindex[index];

    if (n <= flint_primes_small[NEXTPRIME_PRIMES-1])
    {
        if (n == 7) return 7;
        if (n == 11) return 11;
        if (n == 13) return 13;

        while (((n % 7) == 0) || ((n % 11) == 0) || ((n % 13) == 0))
        {
            n += nextmod30[index];
            index = nextindex[index];
        }
        return n;
    }

    moduli = (mp_limb_t *) flint_malloc(NEXTPRIME_PRIMES * sizeof(mp_limb_t));

    for (i = 3; i < NEXTPRIME_PRIMES; i++)
        moduli[i] = (n % flint_primes_small[i]);

    while (1)
    {
        ulong composite = 0;
        ulong diff, acc, pr;

        diff = nextmod30[index];

        /* First check residues */
        for (i = 3; i < NEXTPRIME_PRIMES; i++)
        {
            composite |= (moduli[i] == 0);
            acc = moduli[i] + diff;
            pr = flint_primes_small[i];
            moduli[i] = acc >= pr ? acc - pr : acc;
        }

        if (composite)
        {
            n += diff;
            index = nextindex[index];
            continue;
        }

        if ((!proved && n_is_probabprime(n)) || (proved && n_is_prime(n)))
        {
            break;
        }
        else
        {
            n += diff;
            index = nextindex[index];
        }
    }

    flint_free(moduli);
    return n;
}
