/*
    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2015 Dana Jacobsen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"

static unsigned int nextmod30[] =
{
   1, 6, 5, 4, 3, 2, 1, 4, 3, 2, 1, 2, 1, 4, 3, 2, 1, 2, 1,
   4, 3, 2, 1, 6, 5, 4, 3, 2, 1, 2
};

static unsigned int nextindex[] =
{
   1, 7, 7, 7, 7, 7, 7, 11, 11, 11, 11, 13, 13, 17, 17, 17, 17, 19, 19,
   23, 23, 23, 23, 29, 29, 29, 29, 29, 29, 1
};

/* first 64 primes used for modular arithmetic */
#define N_MODULUS (UWORD(1) << (FLINT_BITS - 1))
#define N_MOD_TAB 64
static const unsigned short n_modular_primes_tab[N_MOD_TAB] = {
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


static mp_limb_t bsearch_uint(mp_limb_t n, const unsigned int *t, int tlen)
{
  int lo = 0;
  int hi = tlen-1;
  while (lo < hi) {
    int mid = lo + (hi-lo)/2;
    if (t[mid] <= n) lo = mid+1;
    else             hi = mid;
  }
  return t[lo];
}

mp_limb_t n_nextprime(mp_limb_t n, int proved)
{
    ulong i, index;

    /* For tiny inputs, bsearch in small primes table */
    if (n < flint_primes_small[FLINT_NUM_PRIMES_SMALL-1])
        return bsearch_uint(n, flint_primes_small, FLINT_NUM_PRIMES_SMALL);

    if (n >= N_MODULUS && n < N_MODULUS + n_modular_primes_tab[N_MOD_TAB-1])
    {
        for (i = 0; i < N_MOD_TAB; i++)
            if (N_MODULUS + n_modular_primes_tab[i] > n)
                return N_MODULUS + n_modular_primes_tab[i];
    }

    if (n >= UWORD_MAX_PRIME)
    {
        flint_throw(FLINT_ERROR, "Exception (n_nextprime). No larger single-limb prime exists.\n");
    }

    index = n % 30;

    do {
      n += nextmod30[index];
      index = nextindex[index];
    } while (!n_is_prime(n));

    return n;
}
