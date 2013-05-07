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

#undef ulong /* prevent clash with standard library */
#include <stdio.h>
#include <stdlib.h>
#define ulong unsigned long
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


mp_limb_t n_nextprime(mp_limb_t n, int proved)
{
    mp_limb_t * moduli;
    unsigned long i, index;

    if (n < 7) 
    {
        if (n < 2)
            return 2;
        n++;
        n |= 1;
        return n;  
    }

    if (n >= ULONG_MAX_PRIME)
    {
        printf("Exception (n_nextprime). No larger single-limb prime exists.\n");
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
        unsigned long composite = 0;
        unsigned long diff, acc, pr;

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
