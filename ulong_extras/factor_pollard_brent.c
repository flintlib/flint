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

mp_limb_t
n_sqr_and_add_a(mp_limb_t y, mp_limb_t a, mp_limb_t n, mp_limb_t ninv,
              mp_limb_t normbits)
{
    mp_limb_t hi, lo;

    y = n_mulmod_preinv(y, y, n, ninv, normbits);
    add_ssaaaa(hi, lo, UWORD(0), y, UWORD(0), a);

    if (hi == 0)
    {
        y = lo;
        if (y > n)
            y -= n;
    }
    else
        sub_ddmmss(hi, y, hi, lo, 0, n);

    return y;
}

int
n_factor_pollard_brent_single(mp_limb_t *factor, mp_limb_t n, mp_limb_t ninv, 
                              mp_limb_t ai, mp_limb_t xi, mp_limb_t normbits,
                              mp_limb_t max_iters)
{
    mp_limb_t iter, i, k, j, minval, m, one_shift_norm, x, y, a, q, ys, subval;
    int ret;

    if (n < 4)
        return 0;

    /* one shifted by normbits, used for comparisions */
    one_shift_norm = UWORD(1) << normbits;
    q = one_shift_norm;
    (*factor) = one_shift_norm;

    y = xi;
    a = ai;

    m = 100;
    iter = 1;
    do {
        x = y;
        k = 0;

        for (i = 0; i < iter; i++)
            y = n_sqr_and_add_a(y, a, n, ninv, normbits);

        do {
            minval = iter - k;
            if (m < minval)
                minval = m;

            ys = y;

            for (i = 0; i < minval; i++)
            {
                y = n_sqr_and_add_a(y, a, n, ninv, normbits);
                if (x > y)
                    subval = x - y;
                else
                    subval = y - x;
                q = n_mulmod_preinv(q, subval, n, ninv, normbits);
            }

            if (q == 0)
                (*factor) = n;
            else
                (*factor) = n_gcd(q, n);
            
            k += m;
            j = ((*factor) == one_shift_norm);
        } while ((k < iter) && (j));

        if (iter > max_iters)
            break;
        iter *= 2;
    }  while (j);

    if ((*factor) == n)
    {
        do {
            ys = n_sqr_and_add_a(ys, a, n, ninv, normbits);
            if (x > ys)
                subval = x - ys;
            else
                subval = ys - x;

            if (q == 0)
                (*factor) = n;
            else
                (*factor) = n_gcd(q, n);
            
            (*factor) = n_gcd(subval, n);
        } while ((*factor) == one_shift_norm);   /* gcd == 1 */
    }

    ret = 1;

    if ((*factor) == one_shift_norm) /* gcd == 1 */
        ret = 0;
    else if ((*factor) == n) /* gcd == n*/
        ret = 0;

    return ret;
}

int
n_factor_pollard_brent(mp_limb_t *factor, flint_rand_t state, mp_limb_t n_in, 
                        mp_limb_t max_tries, mp_limb_t max_iters)
{
    mp_limb_t normbits, a, x, n, ninv, max;
    int ret;

    ret = 0;
    max = n_in -3;    /* 1 <= a <= n - 3 */

    count_leading_zeros(normbits, n_in);
    n = n_in;
    n <<= normbits; 
    ninv = n_preinvert_limb(n);

    while (max_tries--)
    {
        a = n_randint(state, max);
        a += 1;
        max += 1;   /* 1 <= x <= n - 1 */
        x = n_randint(state, max);
        x += 1;

        a <<= normbits;
        x <<= normbits;

        ret = n_factor_pollard_brent_single(factor, n, ninv, a, x, normbits, max_iters);

        if (ret == 1)
        {
            if (normbits)
                (*factor) >>= normbits;
            return 1; 
        }
    }

    return ret;    
}
