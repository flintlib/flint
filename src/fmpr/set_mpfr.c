/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

void
fmpr_set_mpfr(fmpr_t x, const mpfr_t y)
{
    if (!mpfr_regular_p(y))
    {
        if (mpfr_zero_p(y)) fmpr_zero(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) > 0) fmpr_pos_inf(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) < 0) fmpr_neg_inf(x);
        else fmpr_nan(x);
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(fmpr_manref(x));
        fmpz_set_si(fmpr_expref(x), mpfr_get_z_2exp(z, y));
        _fmpz_demote_val(fmpr_manref(x));
        _fmpr_normalise(fmpr_manref(x), fmpr_expref(x), y->_mpfr_prec, FMPR_RND_DOWN);
    }
}

