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
   makes calls to stage I and stage II (two) */

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
    
static
ulong ecm_primorial_two[15] =
{
   UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030), UWORD(510510), UWORD(9699690),
   UWORD(223092870), UWORD(6469693230), UWORD(200560490130), UWORD(7420738134810),
   UWORD(304250263527210), UWORD(13082761331670030), UWORD(614889782588491410)
};

int
fmpz_factor_ecm_two(fmpz_t f, mp_limb_t curves, mp_limb_t B1, mp_limb_t B2,
                    flint_rand_t state, fmpz_t n)
{
    fmpz_t sig, nm8, a24;
    mp_limb_t P, num, maxD, Psq, mmin, mmax, mdiff, prod;
    int i, j, ret;
    ecm_t ecm_inf;

    fmpz_factor_ecm_init(ecm_inf);
    fmpz_init(sig);
    fmpz_init(nm8);
    fmpz_init(a24);
    fmpz_sub_ui(nm8, n, 8);
    
    ret = 0;

    /************************ STAGE I PRECOMPUTATIONS ************************/

    num = n_prime_pi(B1);   /* number of primes under B1 */

    /* compute list of primes under B1 for stage I */
    const mp_limb_t *prime_array = flint_malloc(num * sizeof(mp_limb_t));
    prime_array = n_primes_arr_readonly(num);   

    /************************ STAGE II PRECOMPUTATIONS ***********************/

    maxD = n_sqrt(n_sqrt(B2));

    /* Selecting primorial */

    j = 1;
    while ((j < 15) && (ecm_primorial_two[j] < maxD))
        j += 1;

    P = ecm_primorial_two[j];

    /* compute GCD_table */

    ecm_inf->GCD_table = flint_malloc(P * sizeof(uint8_t));

    for (j = 1; j < P; j += 1)
    {
        if ((j%2) && n_gcd(j, P) == 1)
            ecm_inf->GCD_table[j] = 1;  
        else
            ecm_inf->GCD_table[j] = 0;
    }   
    
    Psq = P*P;
    mmin = (B1) / (Psq);
    mmax = (B2) / (Psq);
    mdiff = mmax - mmin + 1;

    /* compute prime table */

    ecm_inf->prime_table = flint_malloc(mdiff * sizeof(uint8_t*));

    for (i = 0; i < mdiff; i++)
        ecm_inf->prime_table[i] = flint_malloc(Psq * sizeof(uint8_t));

    for (i = 0; i < mdiff; i++)
    {
        for (j = 1; j < Psq; j += 2)
        {
            ecm_inf->prime_table[i][j] = 0;

            if (ecm_inf->GCD_table[j % P] == 1)
            {
                prod = (i + mmin)*Psq + j;

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

        if(fmpz_factor_ecm_stage_II_two(f, B1, B2, P, n, ecm_inf))
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
