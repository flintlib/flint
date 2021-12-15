/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* This is an implementation of the pollard rho algorithm, with a more efficient
   cycle finding algorithm, as proposed by Richard Brent. Details can be found 
   in the paper https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf, pseudocode 
   is available on page 182 of the same paper */

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "mpn_extras.h"

/* Sets y to (y^2 + a) % n */
void
flint_mpn_sqr_and_add_a(mp_ptr y, mp_ptr a, mp_ptr n, mp_limb_t n_size, mp_ptr ninv, 
              mp_limb_t normbits)
{
    mp_limb_t cy;

    flint_mpn_mulmod_preinvn(y, y, y, n_size, n, ninv, normbits);   /* y^2 mod n */
    cy = mpn_add_n(y, y, a, n_size);

    /* Since carry cannot be greater than 1, if carry
       then simply subtract for modulo (a < n, y < n, a + y < 2n).
       If no carry, then check if y > n before subtracting. */

    if (cy)
        mpn_sub_n(y, y, n, n_size);
    else if (mpn_cmp(y, n, n_size) > 0)
            mpn_sub_n(y, y, n, n_size);
}

int
flint_mpn_factor_pollard_brent_single(mp_ptr factor, mp_ptr n, mp_ptr ninv, mp_ptr a, mp_ptr y,
                     mp_limb_t n_size, mp_limb_t normbits, mp_limb_t max_iters)
{     
    /* n_size >= 2, one limb fmpz_t's are passed on to 
       n_factor_pollard_brent in outer funtion      */

    mp_ptr x, q, ys, subval;
    mp_limb_t iter, i, k, minval, m, one_shift_norm, gcdlimbs;
    int ret, j;

    TMP_INIT;
    TMP_START;

    x      = TMP_ALLOC(n_size * sizeof(mp_limb_t));  /* initial value to evaluate f(x) */
    q      = TMP_ALLOC(n_size * sizeof(mp_limb_t));  /* product of gcd's */
    ys     = TMP_ALLOC(n_size * sizeof(mp_limb_t));
    subval = TMP_ALLOC(n_size * sizeof(mp_limb_t));

    /* one shifted by normbits, used for comparisons */
    one_shift_norm = UWORD(1) << normbits;

    /* set factor and q to one (shifted) */
    mpn_zero(q, n_size);        
    mpn_zero(factor, n_size);
    q[0] = one_shift_norm;
    factor[0] = one_shift_norm;

    m = 100;
    iter = 1;
    do {
        flint_mpn_copyi(x, y, n_size);
        k = 0;

        for (i = 0; i < iter; i++)
            flint_mpn_sqr_and_add_a(y, a, n, n_size, ninv, normbits);

        do {
            minval = iter - k;
            if (m < minval)
                minval = m;

            flint_mpn_copyi(ys, y, n_size);

            for (i = 0; i < minval; i++)
            {
                flint_mpn_sqr_and_add_a(y, a, n, n_size, ninv, normbits);
                if (mpn_cmp(x, y, n_size) > 0)
                    mpn_sub_n(subval, x, y, n_size);
                else
                    mpn_sub_n(subval, y, x, n_size);           
                flint_mpn_mulmod_preinvn(q, q, subval, n_size, n, ninv, normbits);  
            }

            /* if q is 0, then gcd is n (gcd(0, a) = a). Not passing through
               flint_mpn_gcd_full due to input paramete restrictions. */
            if (flint_mpn_zero_p(q, n_size) == 0)
                gcdlimbs = flint_mpn_gcd_full(factor, q, n_size, n, n_size);
            else
            {
                flint_mpn_copyi(factor, n, n_size);
                gcdlimbs = n_size;
            }

            k += m;
            j = ((gcdlimbs == 1) && (factor[0] == one_shift_norm));   /* gcd == 1 */
        } while ((k < iter) && j); 

        if (iter > max_iters)   /* max iterations crossed */
            break;

        iter *= 2;
    } while (j);

    /* if gcd == n, could be possible q has all factors of n, start 
       backtracing. if gcd != 1 after backtracing, then at least one 
       factor has been found (can be n) */

    if ((gcdlimbs == n_size) && (mpn_cmp(factor, n, n_size) == 0))  
    {
        do {
            flint_mpn_sqr_and_add_a(ys, a, n, n_size, ninv, normbits);
            if (mpn_cmp(x, ys, n_size) > 0)
                mpn_sub_n(subval, x, ys, n_size);
            else
                mpn_sub_n(subval, ys, x, n_size); 

            if (flint_mpn_zero_p(q, n_size) == 0)
                gcdlimbs = flint_mpn_gcd_full(factor, q, n_size, n, n_size);
            else
            {
                flint_mpn_copyi(factor, n, n_size);
                gcdlimbs = n_size;
            }

            j = ((gcdlimbs == 1) && (factor[0] == one_shift_norm));
        } while (j);   /* gcd == 1 */
    }
    ret = gcdlimbs;

    /* if gcd == 1 or gcd == n, trivial factor found. return 0 */

    if ((gcdlimbs == 1) && (factor[0] == one_shift_norm)) /* gcd == 1 */
        ret = 0;
    else if ((gcdlimbs == n_size && (mpn_cmp(factor, n, n_size) == 0))) /* gcd == n*/
        ret = 0;

    if (ret)
    {
        /* If in case after shifting, "actual" factor has lesser limbs
           than gcdlimbs, then decrease ret by 1. */
        j = n_sizeinbase(factor[gcdlimbs - 1], 2);
        if (normbits >= j)
            ret -= 1;

        if (normbits)
            mpn_rshift(factor, factor, gcdlimbs, normbits);  
    }

    TMP_END;
    
    return ret;
}

int
fmpz_factor_pollard_brent_single(fmpz_t p_factor, fmpz_t n_in, fmpz_t yi, 
                                 fmpz_t ai, mp_limb_t max_iters)
{
    mp_ptr a, y, n, ninv, temp;
    mp_limb_t n_size, normbits, ans, size, cy;
    mp_limb_t al, yl, val, valinv;
    __mpz_struct *fac, *mptr;
    int ret;

    TMP_INIT;

    if (fmpz_is_even(n_in))
    {
       fmpz_set_ui(p_factor, 2);
       return 1;
    }

    n_size = fmpz_size(n_in);

    if (n_size == 1)
    {
        val = fmpz_get_ui(n_in);
        count_leading_zeros(normbits, val);
        val <<= normbits;
        valinv = n_preinvert_limb(val);

        al = fmpz_get_ui(ai);
        yl = fmpz_get_ui(yi);
        al <<= normbits;
        yl <<= normbits;

        ret = n_factor_pollard_brent_single(&ans, val, valinv, al, yl, normbits, max_iters);
        ans >>= normbits;
        fmpz_set_ui(p_factor, ans);
        return ret;
    }

    mptr = COEFF_TO_PTR(*yi);
    temp = COEFF_TO_PTR(*n_in)->_mp_d;
    count_leading_zeros(normbits, temp[n_size - 1]);

    TMP_START;
    a    = TMP_ALLOC(n_size * sizeof(mp_limb_t));
    y    = TMP_ALLOC(n_size * sizeof(mp_limb_t));
    n    = TMP_ALLOC(n_size * sizeof(mp_limb_t));
    ninv = TMP_ALLOC(n_size * sizeof(mp_limb_t));

    /* copying n_in onto n, and normalizing */

    temp = COEFF_TO_PTR(*n_in)->_mp_d;
    count_leading_zeros(normbits, temp[n_size - 1]);
    if (normbits)
        mpn_lshift(n, temp, n_size, normbits);
    else
        flint_mpn_copyi(n, temp, n_size);

    flint_mpn_preinvn(ninv, n, n_size);

    fac = _fmpz_promote(p_factor);
    mpz_realloc2(fac, n_size * FLINT_BITS);
    fac->_mp_size = n_size;

    mpn_zero(a, n_size);
    mpn_zero(y, n_size);

    if (normbits)
    {        
        if ((!COEFF_IS_MPZ(*yi)))
        {
            y[0] = fmpz_get_ui(yi);
            cy = mpn_lshift(y, y, 1, normbits);
            if (cy)
                y[1] = cy;
        }
        else
        {
            mptr = COEFF_TO_PTR(*yi);
            temp = mptr->_mp_d;
            size = mptr->_mp_size;    
            cy = mpn_lshift(y, temp, size, normbits);

            if (cy)
                y[size] = cy;
        }

        if ((!COEFF_IS_MPZ(*ai)))
        {
            a[0] = fmpz_get_ui(ai);
            cy = mpn_lshift(a, a, 1, normbits);
            if (cy)
                a[1] = cy;
        }
        else
        {
            mptr = COEFF_TO_PTR(*ai);
            temp = mptr->_mp_d;
            size = mptr->_mp_size;
            cy = mpn_lshift(a, temp, size, normbits);
            if (cy)
                a[size] = cy;
        }
    }
    else
    {
        temp = COEFF_TO_PTR(*yi)->_mp_d;
        flint_mpn_copyi(y, temp, n_size);
        temp = COEFF_TO_PTR(*ai)->_mp_d;
        flint_mpn_copyi(a, temp, n_size);
    }

    ret = flint_mpn_factor_pollard_brent_single(fac->_mp_d, n, ninv, a, y, n_size, normbits, max_iters);

    if (ret)
    {
        fac->_mp_size = ret;        /* ret is number of limbs of factor found */
        _fmpz_demote_val(p_factor);    
    }

    TMP_END;
    
    return ret;    
}
