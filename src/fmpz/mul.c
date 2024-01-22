/*
    Copyright (C) 2009, 2020 William Hart
    Copyright (C) 2021, 2024 Albin Ahlbäck
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "mpn_extras.h"
#include "fmpz.h"

__GMP_DECLSPEC extern void * (*__gmp_allocate_func)(size_t);
__GMP_DECLSPEC extern void   (*__gmp_free_func)(void *, size_t);

#define GMP_ALLOCATE_LIMBS(nlimbs) __gmp_allocate_func(sizeof(mp_limb_t) * (nlimbs))
#define GMP_FREE_LIMBS(ptr, nlimbs) __gmp_free_func(ptr, sizeof(mp_limb_t) * (nlimbs))

void
fmpz_mul(fmpz_t res, const fmpz_t xp, const fmpz_t yp)
{
    slong xs = *xp;
    slong ys = *yp;
    mpz_ptr mres;
    mp_ptr resd;
    mp_srcptr xd, yd;
    mp_size_t xsz, ysz;
    mp_size_t neg;
    mp_limb_t upperlimb;

    if ((!COEFF_IS_MPZ(xs) && !COEFF_IS_MPZ(ys)) || xs == 0 || ys == 0)
    {
        ulong t[2];
        smul_ppmm(t[1], t[0], xs, ys);
        fmpz_set_signed_uiui(res, t[1], t[0]);
        return;
    }

    mres = _fmpz_promote(res);

    if (!COEFF_IS_MPZ(xs))
        FLINT_SWAP(slong, xs, ys);

    {
        mpz_srcptr mxp = COEFF_TO_PTR(xs);
        xd = mxp->_mp_d;
        xsz = mxp->_mp_size;
    }

    if (!COEFF_IS_MPZ(ys))
    {
        ysz = ys > 0 ? 1 : -1;
        ys = FLINT_ABS(ys);
        yd = (mp_ptr) &ys;
    }
    else
    {
        mpz_srcptr myp = COEFF_TO_PTR(ys);
        yd = myp->_mp_d;
        ysz = myp->_mp_size;
    }

    neg = xsz ^ ysz;
    xsz = FLINT_ABS(xsz);
    ysz = FLINT_ABS(ysz);

    if (xsz < ysz)
    {
        FLINT_SWAP(mp_size_t, xsz, ysz);
        FLINT_SWAP(mp_srcptr, xd, yd);
    }

    resd = mres->_mp_d;

    /* flint_mpn_mul does not work with aliasing apart from when ysz == 1 */
    if (mres->_mp_alloc < xsz + ysz || resd == xd || (resd == yd && ysz != 1))
        resd = GMP_ALLOCATE_LIMBS(xsz + ysz + 1); /* Allocate one extra */

    upperlimb = flint_mpn_mul(resd, xd, xsz, yd, ysz);
    mres->_mp_size = xsz + ysz - (upperlimb == 0);
    if (neg < 0)
        mres->_mp_size = -mres->_mp_size;

    if (resd != mres->_mp_d)
    {
        if (mres->_mp_alloc)
            GMP_FREE_LIMBS(mres->_mp_d, mres->_mp_alloc);
        mres->_mp_d = resd;
        mres->_mp_alloc = xsz + ysz + 1;
    }
}

void
fmpz_mul_ui(fmpz_t res, const fmpz_t xp, ulong ys)
{
    slong xs = *xp;

    if (!COEFF_IS_MPZ(xs) || ys == 0) /* The latter is a hack */
    {
        ulong rd[2];
        ulong abs_xs = FLINT_ABS(xs);

        umul_ppmm(rd[1], rd[0], abs_xs, ys);

        if (rd[1] == 0)
        {
            if (rd[0] <= COEFF_MAX)
            {
                _fmpz_demote(res);
                *res = (xs >= 0) ? rd[0] : -rd[0];
                return;
            }
            else
            {
                mpz_ptr mres = _fmpz_promote(res);
                mres->_mp_d[0] = rd[0];
                mres->_mp_size = (xs >= 0) ? 1 : -1;
                return;
            }
        }
        else
        {
            mpz_ptr mres = _fmpz_promote(res);
            mres->_mp_d[0] = rd[0];
            mres->_mp_d[1] = rd[1];
            mres->_mp_size = (xs >= 0) ? 2 : -2;
            return;
        }
    }
    else
    {
        mpz_ptr mres = _fmpz_promote(res);
        mpz_srcptr mxp = COEFF_TO_PTR(xs);
        mp_size_t sz = mxp->_mp_size;
        mp_size_t abssz = FLINT_ABS(sz);
        mp_limb_t carry;

        if (mres->_mp_alloc <= abssz)
            _mpz_realloc(mres, abssz + 2); /* Allocate one extra */

        /* flint_mpn_mul_N_1 works with aliasing and so does GMP's mpn_mul_1. */
        if (FLINT_HAVE_MUL_FUNC(abssz, 1))
            carry = flint_mpn_mul_func_tab[abssz][1](mres->_mp_d, mxp->_mp_d, &ys);
        else
            carry = mres->_mp_d[abssz] = mpn_mul_1(mres->_mp_d, mxp->_mp_d, abssz, ys);

        mres->_mp_size = abssz + (carry != 0);
        if (sz < 0)
            mres->_mp_size = -mres->_mp_size;
        return;
    }
}
