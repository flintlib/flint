/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "arf.h"

void
arf_set_mpfr(arf_t x, const mpfr_t y)
{
    if (!mpfr_regular_p(y))
    {
        if (mpfr_zero_p(y))
            arf_zero(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) > 0)
            arf_pos_inf(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) < 0)
            arf_neg_inf(x);
        else
            arf_nan(x);
    }
    else
    {
        mp_size_t n = (y->_mpfr_prec + FLINT_BITS - 1) / FLINT_BITS;
        arf_set_mpn(x, y->_mpfr_d, n, y->_mpfr_sign < 0);
        fmpz_set_si(ARF_EXPREF(x), y->_mpfr_exp);
    }
}

