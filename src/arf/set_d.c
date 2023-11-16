/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"

#if FLINT_BITS != 64
# include <mpfr.h>
#endif

#include "double_extras.h"
#include "arf.h"

void
arf_set_d(arf_t x, double v)
{
#if FLINT_BITS == 64
    mp_limb_t h, sign, exp, frac;
    slong real_exp, real_man;
    union { double uf; mp_limb_t ul; } u;

    u.uf = v;
    h = u.ul;
    sign = h >> 63;
    exp = (h << 1) >> 53;
    frac = (h << 12) >> 12;

    if (exp == 0 && frac == 0)
    {
        arf_zero(x);
        return;
    }
    else if (exp == 0x7ff)
    {
        if (frac == 0)
        {
            if (sign)
                arf_neg_inf(x);
            else
                arf_pos_inf(x);
        }
        else
        {
            arf_nan(x);
        }
        return;
    }

    /* handle subnormals */
    if (exp == 0 && frac != 0)
    {
        int exp2;
        v = frexp(v, &exp2);
        u.uf = v;
        h = u.ul;
        sign = h >> 63;
        exp = (h << 1) >> 53;
        frac = (h << 12) >> 12;
        exp += exp2;
    }

    real_exp = exp - 1023 - 52;

    frac |= (UWORD(1) << 52);
    real_man = sign ? (-frac) : frac;

    arf_set_si_2exp_si(x, real_man, real_exp);
#else
    mpfr_t t;
    mp_limb_t tmp[2];

    t->_mpfr_prec = 53;
    t->_mpfr_sign = 1;
    t->_mpfr_exp = 0;
    t->_mpfr_d = tmp;

    mpfr_set_d(t, v, MPFR_RNDD);

    arf_set_mpfr(x, t);
#endif
}

