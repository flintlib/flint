/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "mpn_extras.h"

int
arf_mul_rnd_any(arf_ptr z, arf_srcptr x, arf_srcptr y,
        slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn;
    slong fix;
    int sgnbit, inexact;

    xn = ARF_XSIZE(x);
    yn = ARF_XSIZE(y);
    sgnbit = (xn ^ yn) & 1;
    xn >>= 1;
    yn >>= 1;

    if (yn > xn)
    {
        FLINT_SWAP(arf_srcptr, x, y);
        FLINT_SWAP(mp_size_t, xn, yn);
    }

    /* Either operand is a special value. */
    if (yn == 0)
    {
        arf_mul_special(z, x, y);
        return 0;
    }
    else
    {
        mp_size_t zn, alloc;
        mp_srcptr xptr, yptr;
        mp_ptr tmp;
        ARF_MUL_TMP_DECL

        ARF_GET_MPN_READONLY(xptr, xn, x);
        ARF_GET_MPN_READONLY(yptr, yn, y);

        alloc = zn = xn + yn;
        ARF_MUL_TMP_ALLOC(tmp, alloc)

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xptr, xn, yptr, yn);

        inexact = _arf_set_round_mpn(z, &fix, tmp, zn, sgnbit, prec, rnd);
        _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), fix);
        ARF_MUL_TMP_FREE(tmp, alloc)

        return inexact;
    }
}

