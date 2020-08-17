/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qsieve.h"

/* Array of possible Knuth-Schroeppel multipliers */
static const mp_limb_t multipliers[] = {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15,
                                      17, 19, 21, 22, 23, 26, 29, 30, 31,
                                      33, 34, 35, 37, 38, 41, 42, 43, 47};

/* Number of possible Knuth-Schroeppel multipliers */
#define KS_MULTIPLIERS (sizeof(multipliers)/sizeof(mp_limb_t))

/*
   Try to compute a multiplier k such that there are a lot of small primes
   which are quadratic residues modulo kn. If a small weight of n is found
   during this process it is returned.
*/
mp_limb_t qsieve_knuth_schroeppel(qs_t qs_inf)
{
    float weights[KS_MULTIPLIERS]; /* array of Knuth-Schroeppel weights */
    float best_weight = -10.0f; /* best weight so far */

    ulong i, num_primes, max;
    float logpdivp;
    mp_limb_t nmod8, mod8, p, nmod, pinv, mult;
    int kron, jac;
    n_primes_t iter;
    
    if (fmpz_is_even(qs_inf->n)) /* check 2 is not a factor */
        return 2;

    /* initialise weights for each multiplier k depending on kn mod 8 */
    nmod8 = fmpz_fdiv_ui(qs_inf->n, 8); /* n modulo 8 */

    for (i = 0; i < KS_MULTIPLIERS; i++)
    {
       mod8 = ((nmod8*multipliers[i]) % 8); /* kn modulo 8 */
       weights[i] = 0.34657359; /* ln2/2 */
       if (mod8 == 1) weights[i] *= 4.0;
       if (mod8 == 5) weights[i] *= 2.0;
       weights[i] -= (log((float) multipliers[i]) / 2.0);
    }

    /*
        maximum number of primes to try
        may not exceed number of factor base primes (recall k and 2 and -1 are factor base primes)
    */
    max = FLINT_MIN(qs_inf->ks_primes, qs_inf->num_primes - 3);

    n_primes_init(iter);
    n_primes_next(iter);
    p = n_primes_next(iter);

    for (num_primes = 0; num_primes < max; num_primes++)
    {
        pinv = n_preinvert_limb(p); /* compute precomputed inverse */

        logpdivp = log((float) p) / (float) p; /* log p / p */

        nmod = fmpz_fdiv_ui(qs_inf->n, p);

        if (nmod == 0) return p; /* we found a small factor */

        kron = 1; /* n mod p is even, not handled by n_jacobi */
        while (nmod % 2 == 0)
        {
            if (p % 8 == 3 || p % 8 == 5) kron *= -1;
            nmod /= 2;
        }

        kron *= n_jacobi(nmod, p);
        for (i = 0; i < KS_MULTIPLIERS; i++)
        {
            mult = multipliers[i];
            if (mult >= p)
                mult = n_mod2_preinv(mult, p, pinv); /* k mod p */

            if (mult == 0) weights[i] += logpdivp; /* kn == 0 mod p */
            else
            {
                jac = 1;
                while (mult % 2 == 0) /* k mod p is even, not handled by n_jacobi */
                {
                    if (p % 8 == 3 || p % 8 == 5) jac *= -1;
                    mult /= 2;
                }

                if (kron*jac*n_jacobi(mult, p) == 1) /* kn is a square mod p */
                   weights[i] += 2.0*logpdivp;
            }
        }

        p = n_primes_next(iter);
    }

    n_primes_clear(iter);

    /* search for the multiplier with the best weight and set qs_inf->k */
    for (i = 0; i < KS_MULTIPLIERS; i++)
    {
        if (weights[i] > best_weight)
        {
            best_weight = weights[i];
            qs_inf->k = multipliers[i];
        }
    }

#if QS_DEBUG
    flint_printf("Using Knuth-Schroeppel multiplier %wd\n", qs_inf->k);
#endif

    return 0; /* we didn't find any small factors */
}

