/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz_vec.h"

#define SEARCH_LENGTH 64

void
_fmpz_vec_content(fmpz_t res, const fmpz * vec, slong len)
{
    slong kx;
    ulong t1, t2;
    __mpz_struct * mp;

    /* Strip away all the zeros at the head */
    while (len > 0 && fmpz_is_zero(vec + 0))
    {
        len--;
        vec++;
    }

    /* Strip away all the zeros at the tail */
    while (len > 1 && fmpz_is_zero(vec + len - 1))
        len--;

    if (len <= 3)
    {
        if (len == 0)
            fmpz_zero(res);
        else if (len == 1)
            fmpz_abs(res, vec + 0);
        else if (len == 2)
            fmpz_gcd(res, vec + 0, vec + 1);
        else
            fmpz_gcd3(res, vec + 0, vec + 1, vec + 2);
        return;
    }

    /* The end tends to be +/- 1 */
    if (fmpz_is_pm1(vec + len - 1))
    {
        fmpz_one(res);
        return;
    }

    t2 = FLINT_MIN(len, SEARCH_LENGTH);

    /* Try to get the first small element */
    for (kx = 0; kx < t2; kx++)
    {
        t1 = vec[kx];
        if (t1 == 1 || t1 == -1)
        {
            fmpz_one(res);
            return;
        }
        else if (t1 != 0 && !COEFF_IS_MPZ(t1))
            goto smallfmpz;
    }
    goto unsuresize;

smallfmpz: /* We are sure that the result fits inside a small fmpz */
    t1 = FLINT_ABS(t1); /* t1 is a non-zero small fmpz */
    len -= kx;

    for (; kx > 0; kx--, vec++)             /* Check the zero/mpz elements */
    {
        t2 = *vec;
        if (t2 == 0)
            continue;
        else if (t2 == 1 || t2 == -1 || t1 == 1)
        {
            fmpz_one(res);
            return;
        }

        mp = COEFF_TO_PTR(t2);
        t1 = mpn_gcd_1(mp->_mp_d, FLINT_ABS(mp->_mp_size), t1);
    }
    goto smallfmpz_continuation;

smallfmpz_intermediate: /* If looping in 'unsuresize' detected a small fmpz */
    t2 = FLINT_ABS(t2);
    if (t2 == 1)
    {
        _fmpz_clear_mpz(kx);
        *res = 1;
        return;
    }
    else if (t2 != 0)
    {
        t1 = mpn_gcd_1(mp->_mp_d, mp->_mp_size, t2);
        _fmpz_clear_mpz(kx); /* Clear res/mp, it is no longer needed. */
    }

smallfmpz_continuation:
    for (len--, vec++; len > 0; len--, vec++)   /* And now check the rest */
    {
        t2 = *vec;
        if (t2 == 1 || t2 == -1 || t1 == 1)
        {
            fmpz_one(res);
            return;
        }

        if (!COEFF_IS_MPZ(t2))
        {
            if (t2 == 0)
                continue;
            t2 = FLINT_ABS(t2);
            t1 = mpn_gcd_1(&t2, 1, t1);
        }
        else
        {
            mp = COEFF_TO_PTR(t2);
            t1 = mpn_gcd_1(mp->_mp_d, FLINT_ABS(mp->_mp_size), t1);
        }
    }

    *res = t1;
    return;

unsuresize: /* We are not sure about the result's size */
    if (COEFF_IS_MPZ(*res))    /* If possible, use res's mpz */
    {
        mp = COEFF_TO_PTR(*res);
        kx = *res;  /* kx will be the fmpz-type pointer to mp */
    }
    else
    {
        mp = _fmpz_new_mpz();
        kx = PTR_TO_COEFF(mp);
    }

    t1 = vec[0];
    t2 = vec[1];

    /* At this stage t1 is an mpz and t2 is either an mpz or zero. */
    if (t2 == 0)
        mpz_abs(mp, COEFF_TO_PTR(t1));
    else
        mpz_gcd(mp, COEFF_TO_PTR(t2), COEFF_TO_PTR(t1));

    for (len -= 2, vec += 2; len > 0; len--, vec++)
    {
        t2 = *vec;
        if (t2 == 0)
            continue;
        else if (mp->_mp_size == 1)
        {
            goto ui;
        }
        else if (!COEFF_IS_MPZ(t2))
            goto smallfmpz_intermediate;

        mpz_gcd(mp, COEFF_TO_PTR(t2), mp);
    }

    /* If mp can be represented as a small fmpz */
    if (mp->_mp_size == 1 && mp->_mp_d[0] <= COEFF_MAX)
    {
        *res = mp->_mp_d[0];
        _fmpz_clear_mpz(kx);
        return;
    }

    /* Now mp cannot be represented as a small fmpz. Whether or not res points
     * to mp, it is okay to set res = mp. (If res was an mpz, then res == mp.
     * If not, we simply set res = mp, end of story.) */
    *res = kx;
    return;

ui:
    t1 = mp->_mp_d[0];
    for (; len > 0; len--, vec++)
    {
        t2 = *vec;

        if (t2 == 0)
            continue;
        else if (t1 == 1 || t2 == 1 || t2 == -1)
        {
            _fmpz_clear_mpz(kx);
            *res = 1;
            return;
        }
        else if (!COEFF_IS_MPZ(t2))
        {
            mp->_mp_d[0] = t1;
            mp->_mp_size = 1;
            goto smallfmpz_intermediate;
        }

        mp = COEFF_TO_PTR(t2);
        t1 = mpn_gcd_1(mp->_mp_d, FLINT_ABS(mp->_mp_size), t1);
    }

    if (t1 <= COEFF_MAX)
    {
        _fmpz_clear_mpz(kx);
        *res = t1;
    }
    else
    {
        /* NOTE: After the loop, mp no longer points towards the same object as
         * after the first lines under 'unsuresize'. We reclaim this object by
         * using kx. */
        *res = kx;
        mp = COEFF_TO_PTR(kx);
        mp->_mp_d[0] = t1;
        mp->_mp_size = 1;
    }

    return;
}

#undef SEARCH_LENGTH
