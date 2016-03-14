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
#include "mpn_extras.h"

static
ulong n_ecm_primorial[] =
{
#ifdef FLINT64

    UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030),
    UWORD(510510), UWORD(9699690), UWORD(223092870), UWORD(6469693230), 
    UWORD(200560490130), UWORD(7420738134810), UWORD(304250263527210), 
    UWORD(13082761331670030), UWORD(614889782588491410)
    /* 15 values */

#else

    UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030),
    UWORD(510510), UWORD(9699690)
    /* 9 values */

#endif
};

#ifdef FLINT64
#define num_n_ecm_primorials 15
#else
#define num_n_ecm_primorials 9
#endif

int
fmpz_factor_ecm(fmpz_t f, mp_limb_t curves, mp_limb_t B1, mp_limb_t B2,
                flint_rand_t state, const fmpz_t n_in)
{
    fmpz_t sig, nm8;
    mp_limb_t P, num, maxD, mmin, mmax, mdiff, prod, maxj, n_size, cy;
    int i, j, ret;
    ecm_t ecm_inf;
    __mpz_struct *fac, *mpz_ptr;
    mp_ptr n, mpsig;

    TMP_INIT;

    const mp_limb_t *prime_array;
    n_size = fmpz_size(n_in);

    fmpz_factor_ecm_init(ecm_inf, n_size);

    TMP_START;

    n     = TMP_ALLOC(n_size * sizeof(mp_limb_t));
    mpsig = TMP_ALLOC(n_size * sizeof(mp_limb_t));

    if (n_size == 1)
    {
        ret = n_factor_ecm(&P, curves, B1, B2, state, fmpz_get_ui(n_in));
        fmpz_set_ui(f, P);
        return ret;
    }

    if ((!COEFF_IS_MPZ(* n_in)))
    {
        count_leading_zeros(ecm_inf->normbits, fmpz_get_ui(n_in));
        n[0] = fmpz_get_ui(n_in);
        n[0] <<= ecm_inf->normbits;
    }
    else
    {
        mpz_ptr = COEFF_TO_PTR(* n_in);
        count_leading_zeros(ecm_inf->normbits, mpz_ptr->_mp_d[n_size - 1]);
        mpn_lshift(n, mpz_ptr->_mp_d, n_size, ecm_inf->normbits);
    }

    flint_mpn_preinvn(ecm_inf->ninv, n, n_size);
    ecm_inf->one[0] = UWORD(1) << ecm_inf->normbits;

    fmpz_init(sig);
    fmpz_init(nm8);
    fmpz_sub_ui(nm8, n_in, 8);

    ret = 0;
    fac = _fmpz_promote(f);

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

    ecm_inf->GCD_table = flint_malloc(maxj + 1);

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

        if ((!COEFF_IS_MPZ(*sig)))
        {
            mpsig[0] = fmpz_get_ui(sig);
            cy = mpn_lshift(mpsig, mpsig, 1, ecm_inf->normbits);
            if (cy)
                mpsig[1] = cy;
        }
        else
        {
            mpz_ptr = COEFF_TO_PTR(*sig);

            cy = mpn_lshift(mpsig, mpz_ptr->_mp_d, mpz_ptr->_mp_size, ecm_inf->normbits);
            if (cy)
                mpsig[mpz_ptr->_mp_size] = cy;

        }

        /************************ SELECT CURVE ************************/

        ret = fmpz_factor_ecm_select_curve(fac->_mp_d, mpsig, n, ecm_inf);

        if (ret)
        {
            /* Found factor while selecting curve,
               very very lucky :) */

            mpn_rshift(fac->_mp_d, fac->_mp_d, ret, ecm_inf->normbits);
            MPN_NORM(fac->_mp_d, ret);

            fac->_mp_size = ret;
            _fmpz_demote_val(f);    

            ret = -1;
            goto cleanup;
        }

        if (ret != -1)
        {

            /************************** STAGE I ***************************/

            ret = fmpz_factor_ecm_stage_I(fac->_mp_d, prime_array, num, B1, n, ecm_inf);

            if (ret)
            {
                /* Found factor after stage I */
                mpn_rshift(fac->_mp_d, fac->_mp_d, ret, ecm_inf->normbits);
                MPN_NORM(fac->_mp_d, ret);

                fac->_mp_size = ret;
                _fmpz_demote_val(f);    

                ret = 1;
                goto cleanup;
            }  
            /************************** STAGE II ***************************/

            ret = fmpz_factor_ecm_stage_II(fac->_mp_d, B1, B2, P, n, ecm_inf);

            if (ret)
            {
                /* Found factor after stage I */
                mpn_rshift(fac->_mp_d, fac->_mp_d, ret, ecm_inf->normbits);
                MPN_NORM(fac->_mp_d, ret);
                fac->_mp_size = ret;
                _fmpz_demote_val(f);  

                ret = 2;
                goto cleanup;
            }
        }
    }

    cleanup:


    flint_free(ecm_inf->GCD_table);
    for (i = 0; i < mdiff; i++)
        flint_free(ecm_inf->prime_table[i]);
    flint_free(ecm_inf->prime_table);

    fmpz_factor_ecm_clear(ecm_inf);
    
    TMP_END;

    return ret;
}
