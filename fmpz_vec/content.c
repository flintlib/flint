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

#define SEARCH_LENGTH 24

void
_fmpz_vec_content(fmpz_t res, const fmpz * vec, slong len)
{
    slong kx;
    ulong t1, t2;
    __mpz_struct * mp;

    while (len > 0 && fmpz_is_zero(vec + 0))
    {
        len--;
        vec++;
    }

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
        goto returnone;

    t2 = FLINT_MIN(len, SEARCH_LENGTH);

    /* Try to get the first small element */
    for (kx = 0; kx < t2; kx++)
    {
        t1 = vec[kx];
        if (t1 == 1 || t1 == -1)
            goto returnone;
        else if (t1 != 0 && !COEFF_IS_MPZ(t1))
            goto L1;
    }
    goto L2;

L1: /* We are sure that the result fits inside a small fmpz */
    t1 = FLINT_ABS(vec[kx]);
    len -= kx;

    for (; kx > 0; kx--, vec++)             /* Check the zero/mpz elements */
    {
        t2 = *vec;
        if (t2 == 0)
            continue;
        else if (t2 == 1 || t2 == -1 || t1 == 1)
            goto returnone;

        mp = COEFF_TO_PTR(t2);
        t1 = mpn_gcd_1(mp->_mp_d, FLINT_ABS(mp->_mp_size), t1);
    }
    for (len--, vec++; len > 0; len--, vec++)   /* And now check the rest */
    {
        t2 = *vec;
        if (t2 == 1 || t2 == -1 || t1 == 1)
            goto returnone;

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

L2: /* We are not sure about the results size */
    if (!COEFF_IS_MPZ(*res))    /* If possible, use res's mpz */
        mp = _fmpz_new_mpz();
    else
        mp = COEFF_TO_PTR(*res);

    t1 = vec[0];
    t2 = vec[1];

    if (t1 == 0)
        mpz_set(mp, COEFF_TO_PTR(t2));
    else if (t2 == 0)
        mpz_set(mp, COEFF_TO_PTR(t1));
    else
        mpz_gcd(mp, COEFF_TO_PTR(t2), COEFF_TO_PTR(t1));

    len -= 2;
    vec += 2;

    for (; len > 0; len--, vec++)
    {
        t2 = *vec;
        if (mp->_mp_size == 1)
            goto small;
        else if (!COEFF_IS_MPZ(t2))
            goto presmall;

        mpz_gcd(mp, COEFF_TO_PTR(t2), mp);
    }

    /* _fmpz_demote_val but only the neccessary parts for this function */
    if (mp->_mp_size == 1)
    {
        t1 = mp->_mp_d[0];
        if (t1 <= COEFF_MAX)
        {
            _fmpz_clear_mpz(PTR_TO_COEFF(mp));
            *res = t1;
        }
    }
    else if (mp->_mp_size == 0)
    {
        _fmpz_clear_mpz(PTR_TO_COEFF(mp));
        *res = 0;
    }
    else
        *res = PTR_TO_COEFF(mp);

    return;

presmall:
    if (t2 == 1 || t2 == -1)
        goto returnone;
    else if (t2 != 0)
        mp->_mp_d[0] = mpn_gcd_1(mp->_mp_d, FLINT_ABS(mp->_mp_size), t2);

    len--, vec++;

small:
    t1 = mp->_mp_d[0];
    for (; len > 0; len--, vec++)
    {
        t2 = *vec;

        if (t1 == 1 || t2 == 1 || t2 == -1)
            goto returnone;

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

    fmpz_set_ui(res, t1);
    return;

returnone:
    fmpz_one(res);
    return;
}
