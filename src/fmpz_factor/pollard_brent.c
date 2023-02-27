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

int
fmpz_factor_pollard_brent(fmpz_t p_factor, flint_rand_t state, fmpz_t n_in, 
                          mp_limb_t max_tries, mp_limb_t max_iters)
{
    fmpz_t fa, fy, maxa, maxy;
    mp_ptr a, y, n, ninv, temp;
    mp_limb_t n_size, normbits, ans, val, size, cy;
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
        ret = n_factor_pollard_brent(&ans, state, val, max_tries, max_iters);
        fmpz_set_ui(p_factor, ans);
        return ret;
    }

    fmpz_init2(fa, n_size);
    fmpz_init2(fy, n_size);
    fmpz_init2(maxa, n_size);
    fmpz_init2(maxy, n_size);
    fmpz_sub_ui(maxa, n_in, 3);     /* 1 <= a <= n - 3 */
    fmpz_sub_ui(maxy, n_in, 1);     /* 1 <= y <= n - 1 */

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

    while (max_tries--)
    {
        fmpz_randm(fa, state, maxa);  
        fmpz_add_ui(fa, fa, 1);
        fmpz_randm(fy, state, maxy);
        fmpz_add_ui(fy, fy, 1);

        mpn_zero(a, n_size);
        mpn_zero(y, n_size);

        if (normbits)
        {
            if ((!COEFF_IS_MPZ(*fy)))
            {
                y[0] = fmpz_get_ui(fy);
		cy = mpn_lshift(y, y, 1, normbits);
                if (cy)
                   y[1] = cy;
            }
            else
            {
                mptr = COEFF_TO_PTR(*fy);
                temp = mptr->_mp_d;
                size = mptr->_mp_size;
		cy = mpn_lshift(y, temp, size, normbits);
                if (cy)
                    y[size] = cy;
            }

            if ((!COEFF_IS_MPZ(*fa)))
            {
                a[0] = fmpz_get_ui(fa);
		cy = mpn_lshift(a, a, 1, normbits);
                if (cy)
                    a[1] = cy;
            }
            else
            {
                mptr = COEFF_TO_PTR(*fa);
                temp = mptr->_mp_d;
                size = mptr->_mp_size;
		cy = mpn_lshift(a, temp, size, normbits);
                if (cy)
                    a[size] = cy;
            }
        }
        else
        {
            temp = COEFF_TO_PTR(*fy)->_mp_d;
            flint_mpn_copyi(y, temp, n_size);
            temp = COEFF_TO_PTR(*fa)->_mp_d;
            flint_mpn_copyi(a, temp, n_size);
        }

        ret = flint_mpn_factor_pollard_brent_single(fac->_mp_d, n, ninv, a, y, n_size, normbits, max_iters);

        if (ret)
        {
            fac->_mp_size = ret;        /* ret is number of limbs of factor found */
            _fmpz_demote_val(p_factor);    
            break; 
        }
    }

    fmpz_clear(fa);
    fmpz_clear(fy);
    fmpz_clear(maxa);
    fmpz_clear(maxy);

    TMP_END;
    
    return ret;    
}
