/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2021 Albin Ahlbäck
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "fmpz.h"
#include "mpn_extras.h"

/* This can only be called from fmpz_mul, and assumes
   x and y are not small. */
static void
flint_mpz_mul(mpz_ptr z, mpz_srcptr x, mpz_srcptr y)
{
    mp_size_t xn, yn, zn, sgn;
    mp_srcptr xd, yd;
    mp_ptr zd;
    mp_limb_t top;
    TMP_INIT;

    xn = x->_mp_size;
    yn = y->_mp_size;
    sgn = xn ^ yn;

    xn = FLINT_ABS(xn);
    yn = FLINT_ABS(yn);

    if (xn < yn)
    {
        mpz_srcptr t;
        mp_size_t tn;

        t = x;
        x = y;
        y = t;

        tn = xn;
        xn = yn;
        yn = tn;
    }

    zn = xn + yn;
    if (z->_mp_alloc < zn)
        _mpz_realloc(z, zn);
    zd = z->_mp_d;
    /* Important: read after possibly resizing z, so that the
       pointers are valid in case of aliasing. */
    xd = x->_mp_d;
    yd = y->_mp_d;

    if (xn == yn)
    {
        if (xn == 2)
        {
            mp_limb_t r3, r2, r1, r0;
            FLINT_MPN_MUL_2X2(r3, r2, r1, r0, xd[1], xd[0], yd[1], yd[0]);
            zd[0] = r0;
            zd[1] = r1;
            zd[2] = r2;
            zd[3] = r3;
            zn -= (r3 == 0);
            z->_mp_size = (sgn >= 0) ? zn : -zn;
            return;
        }

        if (xn == 1)
        {
            mp_limb_t hi, lo;
            umul_ppmm(hi, lo,  xd[0], yd[0]);
            zd[0] = lo;
            zd[1] = hi;
            /* The result cannot be 1 limb, because that would
               require a coefficient smaller than COEFF_MAX. */
            FLINT_ASSERT(hi != 0);
            z->_mp_size = (sgn >= 0) ? 2 : -2;
            return;
        }
    }

    /* Unlikely case since operands up to FLINT_BITS-2 bits get
       caught in the fmpz fast path, but we should still optimize
       for this, especially since we don't need to handle aliasing. */
    if (yn == 1)
    {
        if (xn == 2)
        {
            mp_limb_t r2, r1, r0;
            FLINT_MPN_MUL_2X1(r2, r1, r0, xd[1], xd[0], yd[0]);
            zd[0] = r0;
            zd[1] = r1;
            zd[2] = top = r2;
        }
        else
        {
            top = zd[xn] = mpn_mul_1(zd, xd, xn, yd[0]);
        }

        zn -= (top == 0);
        z->_mp_size = (sgn >= 0) ? zn : -zn;
        return;
    }

    TMP_START;

    /* In case of aliasing, we need to copy the input so that
       we do not overwrite it during the multiplication. */
    if (zd == xd)
    {
        mp_ptr tmp = TMP_ALLOC(xn * sizeof(mp_limb_t));
        flint_mpn_copyi(tmp, xd, xn);
        xd = tmp;
    }
    else if (zd == yd)
    {
        mp_ptr tmp = TMP_ALLOC(yn * sizeof(mp_limb_t));
        flint_mpn_copyi(tmp, yd, yn);
        yd = tmp;
    }

    if (x == y)
    {
        flint_mpn_sqr(zd, xd, xn);
        top = zd[zn - 1];
    }
    else
    {
        top = flint_mpn_mul(zd, xd, xn, yd, yn);
    }

    zn -= (top == 0);
    z->_mp_size = (sgn >= 0) ? zn : -zn;

    TMP_END;
}

void
fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    __mpz_struct * mf;
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (!COEFF_IS_MPZ(c1))
    {
        if (!COEFF_IS_MPZ(c2))
        {
            ulong th, tl;
            smul_ppmm(th, tl, c1, c2);
            fmpz_set_signed_uiui(f, th, tl);
            return;
        }
        else if (c1 != 0)
        {
            mf = _fmpz_promote(f);
            flint_mpz_mul_si(mf, COEFF_TO_PTR(c2), c1);
            return;
        }
    }

    if (!COEFF_IS_MPZ(*f))
    {
        if (c1 == 0 || c2 == 0)
        {
            *f = 0;
            return;
        }
        mf = _fmpz_new_mpz();
        (*f) = PTR_TO_COEFF(mf);
    }
    else
    {
        if (c1 == 0 || c2 == 0)
        {
            _fmpz_clear_mpz(*f);
            *f = 0;
            return;
        }

        mf = COEFF_TO_PTR(*f);
    }

    if (!COEFF_IS_MPZ(c2))
        flint_mpz_mul_si(mf, COEFF_TO_PTR(c1), c2);
    else
        flint_mpz_mul(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
}
