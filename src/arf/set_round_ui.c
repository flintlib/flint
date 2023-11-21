/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

/* Top-aligns the single-limb integer value v and rounds it to prec bits.
   Writes inexact, v, exp. Warning: macro without parentheses. */
#define ARF_NORMALISE_ROUND_LIMB(inexact, exp, v, sgnbit, prec, rnd) \
    do { \
        exp = flint_clz(v); \
        v <<= exp; \
        exp = FLINT_BITS - exp; \
        if (prec >= exp) \
        { \
            inexact = 0; \
        } \
        else \
        { \
            mp_limb_t hi_mask, lo_mask, rndn_mask, __t, __u; \
            hi_mask = LIMB_ONES << (FLINT_BITS - prec); \
            __t = v & hi_mask; \
            inexact = (__t != v); \
            if (inexact && rnd != ARF_RND_DOWN) \
            { \
                if (rnd == ARF_RND_NEAR) \
                { \
                    lo_mask = LIMB_ONES >> prec; \
                    rndn_mask = LIMB_ONE << (FLINT_BITS - prec - 1); \
                    __u = v & lo_mask; \
                    if (__u > rndn_mask || (__u == rndn_mask && \
                            (__t << (prec - 1)))) \
                        __t += (LIMB_ONE << (FLINT_BITS - prec)); \
                } \
                else if (arf_rounds_up(rnd, sgnbit)) \
                { \
                    __t += (LIMB_ONE << (FLINT_BITS - prec)); \
                } \
                if (__t == 0) \
                { \
                    __t = LIMB_TOP; \
                    exp++; \
                } \
            } \
            v = __t; \
        } \
    } while (0)

int
_arf_set_round_ui(arf_t x, ulong v, int sgnbit, slong prec, arf_rnd_t rnd)
{
    _fmpz_demote(ARF_EXPREF(x));
    ARF_DEMOTE(x);

    if (v == 0)
    {
        ARF_EXP(x) = ARF_EXP_ZERO;
        ARF_XSIZE(x) = 0;
        return 0;
    }
    else
    {
        int exp, inexact;
        ARF_NORMALISE_ROUND_LIMB(inexact, exp, v, sgnbit, prec, rnd);
        ARF_EXP(x) = exp;
        ARF_XSIZE(x) = ARF_MAKE_XSIZE(1, sgnbit);
        ARF_NOPTR_D(x)[0] = v;
        return inexact;
    }
}

