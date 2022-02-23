/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz.h"

static void
_fmpz_add_mpn_1(fmpz_t f, const mp_limb_t * glimbs, mp_size_t gsz, mp_limb_t x);

static void
_fmpz_sub_mpn_1(fmpz_t f, const mp_limb_t * glimbs, mp_size_t gsz, mp_limb_t x);

void
fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    __mpz_struct * mf;
    slong g1 = *g;
    slong f1 = *f;

    if (!COEFF_IS_MPZ(g1))  /* g is small */
    {
        mp_size_t sz = 2;
        if (g1 >= 0)
        {
            {   /* add with jump if carry */
                ulong tmp = g1;
                g1 += x;
                if (((ulong) g1) < tmp)
                    goto carry;
            }
            if (((ulong) g1) <= COEFF_MAX)
            {
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = g1;
                return;
            }
nocarry:    sz = 1; /* No carry, but result is not a small fmpz */
carry:      if (COEFF_IS_MPZ(f1))
                mf = COEFF_TO_PTR(f1);
            else
            {
                mf = _fmpz_new_mpz();
                *f = PTR_TO_COEFF(mf);
            }
            mf->_mp_size = sz;
            mf->_mp_d[0] = g1;
            mf->_mp_d[1] = 1;
        }
        else
        {
            g1 += x;
            if (((slong) x) >= 0 && g1 <= COEFF_MAX)
            {
                /* x's top bit not set, g1 < 0, and so we can interpret result
                 * as slong. By earlier checking, it is not larger than
                 * COEFF_MAX, and so it must be a small fmpz since x is
                 * positive. */
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = g1;
                return;
            }
            else
            {
                /* If top bit is set in x, the result is going to be positive
                 * and cannot fit in a small fmpz. */

                /* On the other hand, if x's top bit is not set, we can
                 * interpret the result as a slong as g1 is negative. If g1 + x
                 * is larger than COEFF_MAX, then it must be converted to an
                 * mpz. */
                goto nocarry;
            }
        }
    }
    else
    {
        __mpz_struct * mg = COEFF_TO_PTR(g1);
        mp_size_t gsz = mg->_mp_size;
        mp_limb_t * glimbs = mg->_mp_d;

        if (gsz > 0)
            _fmpz_add_mpn_1(f, glimbs, gsz, x);
        else
            _fmpz_sub_mpn_1(f, glimbs, gsz, x);
    }
}

static void
_fmpz_add_mpn_1(fmpz_t f, const mp_limb_t * glimbs, mp_size_t gsz, mp_limb_t x)
{
    __mpz_struct * mf;
    mp_limb_t * flimbs;
    mp_size_t gabssz = FLINT_ABS(gsz);

    if (COEFF_IS_MPZ(*f))
        mf = COEFF_TO_PTR(*f);
    else
    {
        mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
    }
    flimbs = mf->_mp_d;

    if (mf->_mp_alloc < (gabssz + 1))
    {
        mf->_mp_d = realloc(mf->_mp_d, sizeof(mp_limb_t) * (gabssz + 1));
        mf->_mp_alloc = gabssz + 1;
        /* If f and g are aliased, then we need to change glimbs as well. */
        if (flimbs == glimbs)
            glimbs = mf->_mp_d;
        flimbs = mf->_mp_d;
    }

    flimbs[gabssz] = mpn_add_1(flimbs, glimbs, gabssz, x);
    mf->_mp_size = gabssz + flimbs[gabssz];
    if (gsz < 0)
        mf->_mp_size = -mf->_mp_size;
}

static void
_fmpz_sub_mpn_1(fmpz_t f, const mp_limb_t * glimbs, mp_size_t gsz, mp_limb_t x)
{
    __mpz_struct * mf;
    mp_limb_t * flimbs;
    mp_size_t gabssz = FLINT_ABS(gsz);

    if (gabssz == 1)
    {
        if (x <= glimbs[0])
        {
            x = glimbs[0] - x;
L1:         if (x <= COEFF_MAX)
            {
                if (COEFF_IS_MPZ(*f))
                    _fmpz_clear_mpz(*f);
                *f = (gsz > 0) ? x : -x;
            }
            else
            {
                if (COEFF_IS_MPZ(*f))
                    mf = COEFF_TO_PTR(*f);
                else
                {
                    mf = _fmpz_new_mpz();
                    *f = PTR_TO_COEFF(mf);
                }
                mf->_mp_d[0] = x;
                mf->_mp_size = gsz;
            }
        }
        else
        {
            x -= glimbs[0];
            gsz = -gsz; /* Sign change */
            goto L1; /* Reuse code above */
        }
        return;
    }

    if (COEFF_IS_MPZ(*f))
        mf = COEFF_TO_PTR(*f);
    else
    {
        mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
    }
    flimbs = mf->_mp_d;

    if (gabssz == 2)
    {
        sub_ddmmss(flimbs[1], flimbs[0], glimbs[1], glimbs[0], 0, x);
        if (flimbs[1] != 0)
            mf->_mp_size = gsz;
        else if (flimbs[0] > COEFF_MAX)
            mf->_mp_size = (gsz > 0) ? 1 : -1;
        else
        {
            fmpz tmp = *f;
            *f = (gsz > 0) ? flimbs[0] : -flimbs[0];
            _fmpz_clear_mpz(tmp);
        }
    }
    else
    {
        if (mf->_mp_alloc < gabssz)
        {
            /* If f and g cannot be aliased here, since allocation size is
             * always greater than the size's absolute value */
            mf->_mp_d = flimbs = realloc(mf->_mp_d, sizeof(mp_limb_t) * gabssz);
            mf->_mp_alloc = gabssz;
        }
        mpn_sub_1(flimbs, glimbs, gabssz, x);
        mf->_mp_size = gabssz - (flimbs[gabssz - 1] == 0);
        if (gsz < 0)
            mf->_mp_size = -mf->_mp_size;
    }
}

void
fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    __mpz_struct * mf;
    slong g1 = *g;
    slong f1 = *f;

    if (!COEFF_IS_MPZ(g1))  /* g is small */
    {
        mp_size_t sz = -2;
        if (g1 <= 0)
        {
            {   /* "add" with jump if carry */
                ulong tmp = -g1;
                g1 = x - g1; /* g <-- |x| + |g| */
                if (((ulong) g1) < tmp)
                    goto carry;
            }
            if (((ulong) g1) <= COEFF_MAX)
            {
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = -g1;
                return;
            }
nocarry:    sz = -1; /* No carry, but result is not a small fmpz */
carry:      if (COEFF_IS_MPZ(f1))
                mf = COEFF_TO_PTR(f1);
            else
            {
                mf = _fmpz_new_mpz();
                *f = PTR_TO_COEFF(mf);
            }
            mf->_mp_size = sz;
            mf->_mp_d[0] = g1;
            mf->_mp_d[1] = 1;
        }
        else
        {
            g1 = x - g1;
            if (((slong) x) >= 0 && g1 <= COEFF_MAX)
            {
                /* x's top bit not set, g1 > 0, and so we can interpret the
                 * result as slong. And since it is smaller than COEFF_MAX, it
                 * is a small fmpz. */
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = -g1;
                return;
            }
            else
            {
                /* If top bit is set in x, the result is going to be positive
                 * and cannot fit in a small fmpz. */

                /* On the other hand, if x's top bit is not set, but g1 is
                 * larger than COEFF_MAX, then it must be an mpz. */
                goto nocarry;
            }
        }
    }
    else
    {
        __mpz_struct * mg = COEFF_TO_PTR(g1);
        mp_size_t gsz = mg->_mp_size;
        mp_limb_t * glimbs = mg->_mp_d;

        if (gsz > 0)
            _fmpz_sub_mpn_1(f, glimbs, gsz, x);
        else
            _fmpz_add_mpn_1(f, glimbs, gsz, x);
    }
}
