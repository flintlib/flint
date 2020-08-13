/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

static
ulong n_ecm_primorial[] =
{
#ifdef FLINT64

    UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030),
    UWORD(510510), UWORD(9699690), UWORD(223092870), UWORD(6469693230), 
    UWORD(200560490130), UWORD(7420738134810), UWORD(304250263527210), 
    UWORD(13082761331670030), UWORD(614889782588491410)

#else

    UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030),
    UWORD(510510), UWORD(9699690)

#endif
};

#ifdef FLINT64
#define num_n_ecm_primorials 15
#else
#define num_n_ecm_primorials 9
#endif

int
n_factor_ecm(mp_limb_t *f, mp_limb_t curves, mp_limb_t B1, mp_limb_t B2,
             flint_rand_t state, mp_limb_t n)
{
    mp_limb_t P, num, maxD, mmin, mmax, mdiff, prod, maxj, sig;
    int i, j, ret;
    n_ecm_t n_ecm_inf;

    const mp_limb_t *prime_array;

    count_leading_zeros(n_ecm_inf->normbits, n);
    n <<= n_ecm_inf->normbits;
    n_ecm_inf->ninv = n_preinvert_limb(n);
    n_ecm_inf->one = UWORD(1) << n_ecm_inf->normbits;

    ret = 0;

    /************************ STAGE I PRECOMPUTATIONS ************************/

    num = n_prime_pi(B1);   /* number of primes under B1 */

    /* compute list of primes under B1 for stage I */
    prime_array = n_primes_arr_readonly(num);   

    /************************ STAGE II PRECOMPUTATIONS ***********************/

    maxD = n_sqrt(B2);

    /* Selecting primorial */

    j = 1;
    while ((j < num_n_ecm_primorials) && (n_ecm_primorial[j] < maxD))
        j += 1;

    P = n_ecm_primorial[j - 1]; 
    
    mmin = (B1 + (P/2)) / P;
    mmax = ((B2 - P/2) + P - 1)/P;      /* ceil */
    maxj = (P + 1)/2; 
    mdiff = mmax - mmin + 1;

    /* compute GCD_table */

    n_ecm_inf->GCD_table = flint_malloc(maxj + 1);

    for (j = 1; j <= maxj; j += 2)
    {
        if ((j%2) && n_gcd(j, P) == 1)
            n_ecm_inf->GCD_table[j] = 1;  
        else
            n_ecm_inf->GCD_table[j] = 0;
    }  

    /* compute prime table */

    n_ecm_inf->prime_table = flint_malloc(mdiff * sizeof(unsigned char*));

    for (i = 0; i < mdiff; i++)
        n_ecm_inf->prime_table[i] = flint_malloc((maxj + 1) * sizeof(unsigned char));

    for (i = 0; i < mdiff; i++)
    {
        for (j = 1; j <= maxj; j += 2)
        {
            n_ecm_inf->prime_table[i][j] = 0;

            /* if (i + mmin)*D + j
               is prime, mark 1. Can be possibly prime
               only if gcd(j, D) = 1 */

            if (n_ecm_inf->GCD_table[j] == 1)
            {
                prod = (i + mmin)*P + j;
                if (n_is_prime(prod))
                    n_ecm_inf->prime_table[i][j] = 1;

                prod = (i + mmin)*P - j;
                if (n_is_prime(prod))
                    n_ecm_inf->prime_table[i][j] = 1;
            }
        }
    }

    /****************************** TRY "CURVES" *****************************/

    for (j = 0; j < curves; j++)
    {
        sig = n_randint(state, n >> n_ecm_inf->normbits);
        sig = n_addmod(sig, 7, n >> n_ecm_inf->normbits);
        sig <<= n_ecm_inf->normbits;

        /************************ SELECT CURVE ************************/

        if (n_factor_ecm_select_curve(f, sig, n, n_ecm_inf)) 
        {
            /* Found factor while selecting curve,
               very very lucky :) */
            (*f) >>= n_ecm_inf->normbits;
            ret = -1;
            goto cleanup;
        }

        /************************** STAGE I ***************************/

        ret = n_factor_ecm_stage_I(f, prime_array, num, B1, n, n_ecm_inf);

        if (ret)
        {
            /* Found factor after stage I */
            (*f) >>= n_ecm_inf->normbits;
            ret = 1;
            goto cleanup;
        }  

        /************************** STAGE II ***************************/

        ret = n_factor_ecm_stage_II(f, B1, B2, P, n, n_ecm_inf);

        if (ret)
        {
            /* Found factor after stage II */
            (*f) >>= n_ecm_inf->normbits;
            ret = 2;
            goto cleanup;
        }   
    }

    cleanup:

    flint_free(n_ecm_inf->GCD_table);

    for (i = 0; i < mdiff; i++)
        flint_free(n_ecm_inf->prime_table[i]);

    flint_free(n_ecm_inf->prime_table);

    return ret;
}
