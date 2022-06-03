/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2022 Albin Ahlbäck
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz.h"

#define VARIANT 1

void
fmpz_add_ui(fmpz_t res, const fmpz_t x, ulong y)
{
    if (!COEFF_IS_MPZ(*x))
    {
        slong a = *x;

        if (y <= COEFF_MAX)   /* likely case */
        {
#if VARIANT
            slong b = a + y;

            if (b <= COEFF_MAX)
            {
                if (COEFF_IS_MPZ(*res))
                    _fmpz_clear_mpz(*res);
                *res = b;
            }
            else
            {
                __mpz_struct * mpz_res = _fmpz_promote(res);
                flint_mpz_set_si(mpz_res, b);
            }
#else
            fmpz_set_si(res, a + (slong) y);
#endif
        }
        else
        {
            if (a < 0)  /* opposite signs, and y > |a|, so no borrow */
                fmpz_set_ui(res, y + a);
            else
                fmpz_set_uiui(res, ((ulong) a + y) < y, (ulong) a + y);
        }
    }
    else
    {
        mpz_ptr rp;
        mpz_srcptr xp;
        mp_ptr rd;
        mp_srcptr xd;
        mp_size_t xn_signed, xn;
        mp_limb_t cy;

        xp = COEFF_TO_PTR(*x);
        xn_signed = xp->_mp_size;
        xn = FLINT_ABS(xn_signed);

        if (COEFF_IS_MPZ(*res))
            rp = COEFF_TO_PTR(*res);
        else
            rp = _fmpz_promote_val(res);

        if (rp->_mp_alloc < xn + 1)
            _mpz_realloc(rp, xn + 1);

        rd = rp->_mp_d;
        xd = xp->_mp_d;

        if (xn_signed >= 0) /* positive + nonnegative */
        {
            rd[xn] = cy = mpn_add_1(rd, xd, xn, y);
            rp->_mp_size = xn + (cy != 0);
        }
        else if (xn >= 2)  /* negative + nonnegative; result cannot be 0 */
        {
            mpn_sub_1(rd, xd, xn, y);
            xn -= (rd[xn - 1] == 0);
            rp->_mp_size = -xn;

            if (xn == 1 && rd[0] <= COEFF_MAX) /* possible demotion */
            {
                cy = rd[0];
                _fmpz_clear_mpz(*res);
                *res = -(slong) cy;
            }
        }
        else
        {
            /* 1 limb, different signs;
               possible sign change and possible demotion */
            ulong a = xd[0];

            if (y >= a)
            {
                if (y - a <= COEFF_MAX)
                {
                    _fmpz_clear_mpz(*res);
                    *res = y - a;
                }
                else
                {
                    rd[0] = y - a;
                    rp->_mp_size = 1;
                }
            }
            else
            {
                if (a - y <= COEFF_MAX)
                {
                    _fmpz_clear_mpz(*res);
                    *res = -(slong) (a - y);
                }
                else
                {
                    rd[0] = a - y;
                    rp->_mp_size = -1;
                }
            }
        }
    }
}

void
fmpz_sub_ui(fmpz_t res, const fmpz_t x, ulong y)
{
    if (!COEFF_IS_MPZ(*x))
    {
        slong a = *x;

        if (y <= COEFF_MAX)
        {
#if VARIANT
            slong b = a - y;

            if (b >= COEFF_MIN)
            {
                if (COEFF_IS_MPZ(*res))
                    _fmpz_clear_mpz(*res);
                *res = b;
            }
            else
            {
                __mpz_struct * mpz_res = _fmpz_promote(res);
                flint_mpz_set_si(mpz_res, b);
            }
#else
            fmpz_set_si(res, a - (slong) y);
#endif
        }
        else
        {
            if (a > 0)  /* opposite signs, and y > |a|, so no borrow */
                fmpz_neg_ui(res, y - a);
            else
                fmpz_neg_uiui(res, ((ulong) (-a) + y) < y, (ulong) (-a) + y);
        }
    }
    else
    {
        mpz_ptr rp;
        mpz_srcptr xp;
        mp_ptr rd;
        mp_srcptr xd;
        mp_size_t xn_signed, xn;
        mp_limb_t cy;

        xp = COEFF_TO_PTR(*x);
        xn_signed = xp->_mp_size;
        xn = FLINT_ABS(xn_signed);

        if (COEFF_IS_MPZ(*res))
            rp = COEFF_TO_PTR(*res);
        else
            rp = _fmpz_promote_val(res);

        if (rp->_mp_alloc < xn + 1)
            _mpz_realloc(rp, xn + 1);

        rd = rp->_mp_d;
        xd = xp->_mp_d;

        if (xn_signed <= 0) /* positive + nonnegative */
        {
            rd[xn] = cy = mpn_add_1(rd, xd, xn, y);
            rp->_mp_size = -(xn + (cy != 0));
        }
        else if (xn >= 2)  /* negative + nonnegative; result cannot be 0 */
        {
            mpn_sub_1(rd, xd, xn, y);
            xn -= (rd[xn - 1] == 0);
            rp->_mp_size = xn;

            if (xn == 1 && rd[0] <= COEFF_MAX) /* possible demotion */
            {
                cy = rd[0];
                _fmpz_clear_mpz(*res);
                *res = (slong) cy;
            }
        }
        else
        {
            /* 1 limb, different signs;
               possible sign change and possible demotion */
            ulong a = xd[0];

            if (y >= a)
            {
                if (y - a <= COEFF_MAX)
                {
                    _fmpz_clear_mpz(*res);
                    *res = -(slong) (y - a);
                }
                else
                {
                    rd[0] = y - a;
                    rp->_mp_size = -1;
                }
            }
            else
            {
                if (a - y <= COEFF_MAX)
                {
                    _fmpz_clear_mpz(*res);
                    *res = a - y;
                }
                else
                {
                    rd[0] = a - y;
                    rp->_mp_size = 1;
                }
            }
        }
    }
}
