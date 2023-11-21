/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "mpn_extras.h"

int
arf_fma(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn, zn, tn, alloc;
    mp_srcptr xptr, yptr, zptr;
    mp_ptr tptr, tptr2;
    fmpz_t texp;
    slong shift;
    int tsgnbit, inexact;
    ARF_MUL_TMP_DECL

    if (arf_is_special(x) || arf_is_special(y) || arf_is_special(z))
    {
        if (arf_is_zero(z))
        {
            return arf_mul(res, x, y, prec, rnd);
        }
        else if (arf_is_finite(x) && arf_is_finite(y))
        {
            return arf_set_round(res, z, prec, rnd);
        }
        else
        {
            /* todo: speed up */
            arf_t t;
            arf_init(t);
            arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_add(res, z, t, prec, rnd);
            arf_clear(t);
            return inexact;
        }
    }

    tsgnbit = ARF_SGNBIT(x) ^ ARF_SGNBIT(y);
    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    ARF_GET_MPN_READONLY(zptr, zn, z);

    fmpz_init(texp);

    _fmpz_add2_fast(texp, ARF_EXPREF(x), ARF_EXPREF(y), 0);
    shift = _fmpz_sub_small(ARF_EXPREF(z), texp);

    alloc = tn = xn + yn;
    ARF_MUL_TMP_ALLOC(tptr2, alloc)
    tptr = tptr2;

    FLINT_MPN_MUL_WITH_SPECIAL_CASES(tptr, xptr, xn, yptr, yn);

    tn -= (tptr[0] == 0);
    tptr += (tptr[0] == 0);

    if (shift >= 0)
        inexact = _arf_add_mpn(res, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
            tptr, tn, tsgnbit, shift, prec, rnd);
    else
        inexact = _arf_add_mpn(res, tptr, tn, tsgnbit, texp,
            zptr, zn, ARF_SGNBIT(z), -shift, prec, rnd);

    ARF_MUL_TMP_FREE(tptr2, alloc)
    fmpz_clear(texp);

    return inexact;
}

