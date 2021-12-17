/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

#if FLINT64
#define LARGEST_ULONG_PRIMORIAL 52
#else
#define LARGEST_ULONG_PRIMORIAL 28
#endif

/* Only those with odd index */
const ulong ULONG_PRIMORIALS[] = 
{
    UWORD(6),
    UWORD(30),
    UWORD(210),
    UWORD(210),
    UWORD(2310),
    UWORD(30030),
    UWORD(30030),
    UWORD(510510),
    UWORD(9699690),
    UWORD(9699690),
    UWORD(223092870),
    UWORD(223092870),
    UWORD(223092870),
#if FLINT64
    UWORD(6469693230),
    UWORD(200560490130),
    UWORD(200560490130),
    UWORD(200560490130),
    UWORD(7420738134810),
    UWORD(7420738134810),
    UWORD(304250263527210),
    UWORD(13082761331670030),
    UWORD(13082761331670030),
    UWORD(614889782588491410),
    UWORD(614889782588491410),
    UWORD(614889782588491410)
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
        result[0] = UWORD(1);
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

void
fmpz_primorial(fmpz_t res, ulong n)
{
    mp_size_t len, pi;
    ulong bits;
    __mpz_struct * mres;
    const mp_limb_t * primes;

    if (n <= LARGEST_ULONG_PRIMORIAL)
    {
        if (n <= 2)
            fmpz_set_ui(res, 1 + (n==2));
        else
            fmpz_set_ui(res, ULONG_PRIMORIALS[(n-1)/2-1]);
        return;
    }

    pi = n_prime_pi(n);

    primes = n_primes_arr_readonly(pi);
    bits = FLINT_BIT_COUNT(primes[pi - 1]);
    
    mres = _fmpz_promote(res);
    mpz_realloc2(mres, pi*bits);
    
    len = mpn_prod_limbs(mres->_mp_d, primes, pi, bits);
    mres->_mp_size = len;
}

