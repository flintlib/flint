/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arf.h"

#if FLINT_HAVE_MPN_MULHIGH_NORMALISED
/* NOTE: Assumes no special values */
void
arf_mul_rnd_sloppy(arf_ptr z, arf_srcptr x, arf_srcptr y)
{
    mp_size_t xn, yn, zn;
    mp_srcptr xp, yp;
    mp_ptr zp;
    int sgnbit;
    struct mp_limb_pair_t mulret;

    xn = ARF_XSIZE(x);
    yn = ARF_XSIZE(y);
    zn = ARF_XSIZE(z);
    sgnbit = (xn ^ yn) & 1;
    xn >>= 1;
    yn >>= 1;
    zn >>= 1;

    if (yn > xn)
    {
        FLINT_SWAP(arf_srcptr, x, y);
        FLINT_SWAP(mp_size_t, xn, yn);
    }

    if (yn > 2)
    {
        xp = ARF_PTR_D(x);
        yp = ARF_PTR_D(y);
    }
    else
    {
        yp = ARF_NOPTR_D(y);
        if (xn <= 2)
            xp = ARF_NOPTR_D(x);
        else
            xp = ARF_PTR_D(x);
    }

    xp += xn - yn; /* Make xp have the same length as yp */

    if (zn < yn)
    {
        _arf_promote(z, yn);
        zp = ARF_PTR_D(z);
    }
    else if (zn > 2 && yn <= 2)
    {
        _arf_demote(z);
        zp = ARF_NOPTR_D(z);
    }

    if (!FLINT_HAVE_MULHIGH_NORMALISED_FUNC(yn) && (zp == xp || zp == yp))
    {
        /* These cases do not allow aliasing */
        mp_ptr tmp = flint_malloc(sizeof(mp_limb_t) * yn);
        mulret = flint_mpn_mulhigh_normalised(tmp, xp, yp, yn);
        flint_mpn_copyi(zp, tmp, yn);
        flint_free(tmp);
    }
    else
    {
        /* Here we actually allow aliasing. */
        mulret = flint_mpn_mulhigh_normalised(zp, xp, yp, yn);
    }

    ARF_XSIZE(z) = ARF_MAKE_XSIZE(yn, sgnbit);
    _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), mulret.m2);
}
#else
typedef int fileisempty;
#endif
