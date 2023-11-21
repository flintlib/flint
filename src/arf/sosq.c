/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "mpn_extras.h"

int
arf_sosq(arf_t res, const arf_t a, const arf_t b, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(a) || arf_is_special(b))
    {
        if (arf_is_zero(a))
            return arf_mul(res, b, b, prec, rnd);

        if (arf_is_zero(b))
            return arf_mul(res, a, a, prec, rnd);

        if (arf_is_nan(a) || arf_is_nan(b))
            arf_nan(res);
        else
            arf_pos_inf(res);

        return 0;
    }
    else
    {
        mp_srcptr ap, bp;
        int inexact;
        mp_ptr tmp, aap, bbp;
        mp_size_t an, bn, aan, bbn, alloc;
        slong shift;
        fmpz_t texp, uexp;
        ARF_MUL_TMP_DECL

        ARF_GET_MPN_READONLY(ap, an, a);
        ARF_GET_MPN_READONLY(bp, bn, b);

        fmpz_init(texp);
        fmpz_init(uexp);

        _fmpz_add2_fast(texp, ARF_EXPREF(a), ARF_EXPREF(a), 0);
        _fmpz_add2_fast(uexp, ARF_EXPREF(b), ARF_EXPREF(b), 0);
        shift = _fmpz_sub_small(texp, uexp);

        aan = 2 * an;
        bbn = 2 * bn;
        alloc = aan + bbn;

        ARF_MUL_TMP_ALLOC(tmp, alloc)
        aap = tmp;
        bbp = tmp + aan;

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(aap, ap, an, ap, an)
        aan -= (aap[0] == 0);
        aap += (aap[0] == 0);

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(bbp, bp, bn, bp, bn)
        bbn -= (bbp[0] == 0);
        bbp += (bbp[0] == 0);

        if (shift >= 0)
            inexact = _arf_add_mpn(res, aap, aan, 0, texp,
                                        bbp, bbn, 0, shift, prec, rnd);
        else
            inexact = _arf_add_mpn(res, bbp, bbn, 0, uexp,
                                        aap, aan, 0, -shift, prec, rnd);

        ARF_MUL_TMP_FREE(tmp, alloc)
        fmpz_clear(texp);
        fmpz_clear(uexp);

        return inexact;
    }
}

