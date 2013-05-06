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

    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "arith.h"
#include "ulong_extras.h"

#if FLINT64
#define LARGEST_ULONG_PRIMORIAL 52
#else
#define LARGEST_ULONG_PRIMORIAL 28
#endif

/* Only those with odd index */
const ulong ULONG_PRIMORIALS[] = 
{
    6UL,30UL,210UL,210UL,2310UL,30030UL,30030UL,510510UL,9699690UL,9699690UL,
    223092870UL,223092870UL,223092870UL
#if FLINT64
    ,6469693230UL,200560490130UL,200560490130UL,200560490130UL,7420738134810UL,
    7420738134810UL,304250263527210UL,13082761331670030UL,13082761331670030UL,
    614889782588491410UL, 614889782588491410UL, 614889782588491410UL
#endif
};


#define PROD_LIMBS_DIRECT_CUTOFF 50

mp_size_t mpn_prod_limbs_direct(mp_limb_t * result, const mp_limb_t * factors,
    mp_size_t n)
{
    mp_size_t k, len;
    mp_limb_t top;
    if (n < 1)
    {
        result[0] = 1UL;
        return 1;
    }
    result[0] = factors[0];
    len = 1;
    for (k=1; k<n; k++)
    {
        top = mpn_mul_1(result, result, len, factors[k]);
        if (top)
        {
            result[len] = top;
            len++;
        }
    }
    return len;
}

mp_size_t mpn_prod_limbs_balanced(mp_limb_t * result, mp_limb_t * scratch,
                             const mp_limb_t * factors, mp_size_t n, ulong bits)
{
    mp_size_t an, bn, alen, blen, len;
    mp_limb_t top;

    if (n < PROD_LIMBS_DIRECT_CUTOFF)
        return mpn_prod_limbs_direct(result, factors, n);

    an = n/2;
    bn = n - an;
    
    alen = mpn_prod_limbs_balanced(scratch, result, factors, an, bits);
    blen = mpn_prod_limbs_balanced(scratch + alen, result, factors + an, bn, bits);
    len = alen + blen;

    if (alen <= blen)
        top = mpn_mul(result, scratch + alen, blen, scratch, alen);
    else
        top = mpn_mul(result, scratch, alen, scratch + alen, blen);

    if (!top)
        len--;
    
    return len;
}

/*
    Set result to the product of the given factors, return the
    length of the result. It is assumed that no factors are zero.
    bits must be set to some bound on the bit size of the entries
    in factors. If no bound is known, simply use FLINT_BITS.
*/
mp_size_t mpn_prod_limbs(mp_limb_t * result, const mp_limb_t * factors,
    mp_size_t n, ulong bits)
{
    mp_size_t len, limbs;
    mp_limb_t * scratch;
    
    if (n < PROD_LIMBS_DIRECT_CUTOFF)
        return mpn_prod_limbs_direct(result, factors, n);

    limbs = (n * bits - 1)/FLINT_BITS + 2; 

    scratch = flint_malloc(sizeof(mp_limb_t) * limbs);
    len = mpn_prod_limbs_balanced(result, scratch, factors, n, bits);
    flint_free(scratch);
    
    return len;
}

void arith_primorial(fmpz_t res, long n)
{
    mp_size_t len, pi;
    ulong bits;
    __mpz_struct * mpz_ptr;

    if (n <= LARGEST_ULONG_PRIMORIAL)
    {
        if (n <= 2)
            fmpz_set_ui(res, 1 + (n==2));
        else
            fmpz_set_ui(res, ULONG_PRIMORIALS[(n-1)/2-1]);
        return;
    }

    pi = n_prime_pi(n);
    
    n_compute_primes(pi);
    bits = FLINT_BIT_COUNT(flint_primes[pi - 1]);
    
    mpz_ptr = _fmpz_promote(res);
    mpz_realloc2(mpz_ptr, pi*bits);
    
    len = mpn_prod_limbs(mpz_ptr->_mp_d, flint_primes, pi, bits);
    mpz_ptr->_mp_size = len;
}
