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

/* Outer wrapper for ECM 
   makes calls to stage I and stage II (one) */

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
    
static
ulong ecm_primorial_one[15] =
{
   UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030), UWORD(510510), UWORD(9699690),
   UWORD(223092870), UWORD(6469693230), UWORD(200560490130), UWORD(7420738134810),
   UWORD(304250263527210), UWORD(13082761331670030), UWORD(614889782588491410)
};

int
fmpz_factor_ecm_one(fmpz_t f, mp_limb_t curves, mp_limb_t B1, mp_limb_t B2,
                    flint_rand_t state, fmpz_t n)
{
    fmpz_t sig, nm8, a24;
    mp_limb_t P, num, maxD, mmin, mmax, mdiff, prod, maxj;
    int i, j, ret;
    ecm_t ecm_inf;

    const mp_limb_t *prime_array;
    
    fmpz_factor_ecm_init(ecm_inf);
    fmpz_init(sig);
    fmpz_init(nm8);
    fmpz_init(a24);
    fmpz_sub_ui(nm8, n, 8);
    
    ret = 0;

    /************************ STAGE I PRECOMPUTATIONS ************************/

    num = n_prime_pi(B1);   /* number of primes under B1 */

    /* compute list of primes under B1 for stage I */
    prime_array = n_primes_arr_readonly(num);   

    /************************ STAGE II PRECOMPUTATIONS ***********************/

    maxD = n_sqrt(B2);

    /* Selecting primorial */

    j = 1;
    while ((j < 15) && (ecm_primorial_one[j] < maxD))
        j += 1;

    P = ecm_primorial_one[j - 1]; 
    
    mmin = (B1 + (P/2)) / P;
    mmax = ((B2 - P/2) + P - 1)/P;      /* ceil */
    maxj = (P + 1)/2; 
    mdiff = mmax - mmin + 1;

    /* compute GCD_table */

    ecm_inf->GCD_table = flint_malloc((maxj + 1) * sizeof(unsigned char));

    for (j = 1; j <= maxj; j += 2)
    {
        if ((j%2) && n_gcd(j, P) == 1)
            ecm_inf->GCD_table[j] = 1;  
        else
            ecm_inf->GCD_table[j] = 0;
    }  

    /* compute prime table */

    ecm_inf->prime_table = flint_malloc(mdiff * sizeof(unsigned char*));

    for (i = 0; i < mdiff; i++)
        ecm_inf->prime_table[i] = flint_malloc((maxj + 1) * sizeof(unsigned char));

    for (i = 0; i < mdiff; i++)
    {
        for (j = 1; j <= maxj; j += 2)
        {
            ecm_inf->prime_table[i][j] = 0;

            /* if (i + mmin)*D + j
               is prime, mark 1. Can be possibly prime
               only if gcd(j, D) = 1 */

            if (ecm_inf->GCD_table[j] == 1)
            {
                prod = (i + mmin)*P + j;
                if (n_is_prime(prod))
                    ecm_inf->prime_table[i][j] = 1;

                prod = (i + mmin)*P - j;
                if (n_is_prime(prod))
                    ecm_inf->prime_table[i][j] = 1;
            }
        }
    }

    /****************************** TRY "CURVES" *****************************/

    for (j = 0; j < curves; j++)
    {
        fmpz_randm(sig, state, nm8);
        fmpz_add_ui(sig, sig, 7);


        /************************ SELECT CURVE ************************/

        if (fmpz_factor_ecm_select_curve(f, sig, n, ecm_inf)) 
        {
            /* Found factor while selecting curve,
               very very lucky :) */
            ret = 1;
            goto cleanup;
        }        

        /************************** STAGE I ***************************/

        if (fmpz_factor_ecm_stage_I(f, prime_array, num, B1, n, ecm_inf))
        {
            /* Found factor after stage I */
            ret = 1;
            goto cleanup;
        }

        /************************** STAGE II ***************************/

        if(fmpz_factor_ecm_stage_II_one(f, B1, B2, P, n, ecm_inf))
        {
            /* Found factor after stage II */
            ret = 1;
            goto cleanup;
        }        
    }

    cleanup:

    fmpz_clear(sig);
    fmpz_clear(nm8);
    fmpz_clear(a24);
    fmpz_factor_ecm_clear(ecm_inf);

    return ret;
}
