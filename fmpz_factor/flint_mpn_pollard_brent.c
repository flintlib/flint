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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"
    
void
mpn_sqr_and_add_a(mp_ptr y, mp_ptr a, mp_ptr n, mp_limb_t n_size, mp_ptr ninv, mp_limb_t normbits)
{
    mp_ptr temp, tempn;
    mp_limb_t new_size, cy;

    new_size = n_size + 1;

    temp = flint_malloc((new_size) * sizeof(mp_limb_t));
    tempn = flint_malloc((new_size) * sizeof(mp_limb_t));

    mpn_copyi(tempn, n, n_size);
    tempn[n_size] = UWORD(0);

    flint_mpn_mulmod_preinvn(y, y, y, n_size, n, ninv, normbits);
    cy = mpn_add_n(temp, y, a, n_size);

    if (cy)
    {
        temp[n_size] = cy;        
        mpn_sub_n(temp, temp, tempn, new_size);
        mpn_copyi(y, temp, n_size);
    }
    else
    {
        if (mpn_cmp(temp, n, n_size) > 0)
            mpn_sub_n(y, temp, n, n_size);
        else
            mpn_copyi(y, temp, n_size);
    }

    flint_free(temp);
    flint_free(tempn);
}

int
flint_mpn_factor_pollard_brent_single(mp_ptr gcdval, mp_ptr n, mp_ptr ninv, mp_ptr ai, 
                                 mp_ptr xi, mp_limb_t n_size, mp_limb_t normbits, 
                                 mp_limb_t max_iters)
{

    mp_ptr x, y, a, q, ys, subval;
    mp_limb_t iter, i, k, j, minval, m, one_shift_norm, gcdlimbs;
    int ret;

    x      = flint_malloc(n_size * sizeof(mp_limb_t));
    y      = flint_malloc(n_size * sizeof(mp_limb_t));
    a      = flint_malloc(n_size * sizeof(mp_limb_t));
    q      = flint_malloc(n_size * sizeof(mp_limb_t));
    ys     = flint_malloc(n_size * sizeof(mp_limb_t));
    subval = flint_malloc(n_size * sizeof(mp_limb_t));

    /* one shifted by normbits, used for comparisions */
    one_shift_norm = UWORD(1) << normbits;

    mpn_zero(q, n_size);
    mpn_zero(gcdval, n_size);

    q[0] = one_shift_norm;
    gcdval[0] = one_shift_norm;

    mpn_copyi(y, xi, n_size);
    mpn_copyi(a, ai, n_size);

    m = 100;
    iter = 1;
    do {
        mpn_copyi(x, y, n_size);
        k = 0;

        for (i = 0; i < iter; i++)
            mpn_sqr_and_add_a(y, a, n, n_size, ninv, normbits);

        do {
            minval = iter - k;
            if (m < minval)
                minval = m;

            mpn_copyi(ys, y, n_size);

            for (i = 0; i < minval; i++)
            {
                mpn_sqr_and_add_a(y, a, n, n_size, ninv, normbits);
                if (mpn_cmp(x, y, n_size) > 0)
                    mpn_sub_n(subval, x, y, n_size);
                else
                    mpn_sub_n(subval, y, x, n_size);           
                flint_mpn_mulmod_preinvn(q, q, subval, n_size, n, ninv, normbits);  
            }

            if (flint_mpn_zero_p(q, n_size))
            {
                ret = 0;
                goto cleanup;
            }
            gcdlimbs = flint_mpn_gcd_full(gcdval, q, n_size, n, n_size);

            k += m;

            j = ((gcdlimbs == 1) && gcdval[0] == one_shift_norm);

        } while ((k < iter) && j);

        if (iter > max_iters)
            break;

        iter *= 2;
    }  while (j);

    if (gcdlimbs == n_size && !mpn_cmp(gcdval, n, n_size))
    {
        do {
            mpn_sqr_and_add_a(ys, a, n, n_size, ninv, normbits);
            mpn_sub_n(subval, x, ys, n_size);

            if (flint_mpn_zero_p(subval, n_size))
            {
                ret = 0;
                goto cleanup;
            }

            gcdlimbs = flint_mpn_gcd_full(gcdval, subval, n_size, n, n_size);
        } while ((gcdlimbs == 1) && gcdval[0] == one_shift_norm);
    }
    ret = 1;

    if ((gcdlimbs == 1 && gcdval[0] == one_shift_norm))
        ret = 0;
    else if ((gcdlimbs == n_size && !mpn_cmp(gcdval, n, n_size)))
        ret = 0;

    if (ret)
    {
        if (normbits)
        {
            if (n_size == 1)
                gcdval[0] >>= normbits;
            else
                mpn_rshift(gcdval, gcdval, gcdlimbs, normbits);
        }
    }

    cleanup:

    flint_free(x);
    flint_free(y);
    flint_free(a);
    flint_free(q);
    flint_free(ys);
    flint_free(subval);

    return ret;
}

int
flint_mpn_factor_pollard_brent(fmpz_t p_factor, flint_rand_t state, fmpz_t n_in, 
                        mp_limb_t max_tries, mp_limb_t max_iters)
{
    fmpz_t fa, fx, maxa, maxx;
    mp_ptr a, x, n, ninv, temp, tempa, tempx;
    mp_limb_t n_size, normbits;
    int ret;

    n_size = fmpz_size(n_in);

    fmpz_init(fa);
    fmpz_init(fx);
    fmpz_init(maxa);
    fmpz_init(maxx);
    fmpz_init(p_factor);

    ret = 0;

    fmpz_sub_ui(maxa, n_in, 3);
    fmpz_sub_ui(maxx, n_in, 1);


    a = flint_malloc(n_size * sizeof(mp_limb_t));
    x = flint_malloc(n_size * sizeof(mp_limb_t));
    n = flint_malloc(n_size * sizeof(mp_limb_t));
    ninv = flint_malloc(n_size * sizeof(mp_limb_t));

    /* copying n_in onto n, and normalizing */
    if (n_size == 1)
    {
        n[0] = fmpz_get_ui(n_in);
        count_leading_zeros(normbits, n[0]);
        n[0] <<= normbits; 
    }
    else
    {
        temp = COEFF_TO_PTR(*n_in)->_mp_d;
        count_leading_zeros(normbits, temp[n_size - 1]);
        if (normbits)
            mpn_lshift(n, temp, n_size, normbits);
        else
            mpn_copyi(n, temp, n_size);
    }

    flint_mpn_preinvn(ninv, n, n_size);

    __mpz_struct *fac = _fmpz_promote(p_factor);
    mpz_realloc2(fac, n_size * FLINT_BITS);
    fac->_mp_size = n_size;
    while (max_tries--)
    {
        fmpz_randm(fa, state, maxa);  
        fmpz_add_ui(fa, fa, 1);
        fmpz_randm(fx, state, maxx);
        fmpz_add_ui(fx, fx, 1);

        if (n_size == 1)
        {
            x[0] = fmpz_get_ui(fx);
            x[0] <<= normbits; 
            a[0] = fmpz_get_ui(fa);
            a[0] <<= normbits; 
        }
        else
        {
            if (normbits)
            {
                if (fmpz_size(fx) == 1)
                {
                    x[0] = fmpz_get_ui(fx);
                    x[1] = UWORD(0);
                    mpn_lshift(x, x, n_size, normbits);
                }
                else
                {
                    tempx = COEFF_TO_PTR(*fx)->_mp_d;
                    mpn_lshift(x, tempx, n_size, normbits);
                }

                if (fmpz_size(fa) == 1)
                {
                    a[0] = fmpz_get_ui(fa);
                    a[1] = UWORD(0);
                    mpn_lshift(a, a, n_size, normbits);
                }
                else
                {
                    tempa = COEFF_TO_PTR(*fa)->_mp_d;
                    mpn_lshift(a, tempa, n_size, normbits);
                }
            }
            else
            {
                temp = COEFF_TO_PTR(*fx)->_mp_d;
                mpn_copyi(x, temp, n_size);
                temp = COEFF_TO_PTR(*fa)->_mp_d;
                mpn_copyi(a, temp, n_size);
            }
        }

        ret = flint_mpn_factor_pollard_brent_single(p_factor, n, ninv, a, x, n_size, normbits, max_iters);

        if (ret == 1)
            break; 
    }

    fmpz_clear(fa);
    fmpz_clear(fx);
    fmpz_clear(maxa);
    fmpz_clear(maxx);

    flint_free(a);
    flint_free(x);
    flint_free(n);
    flint_free(ninv);

    return ret;    
}
