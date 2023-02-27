/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

slong _fmpr_add_eps(fmpr_t z, const fmpr_t x, int sign, slong prec, fmpr_rnd_t rnd)
{
    slong bits, shift;
    int xsign;

    xsign = fmpz_sgn(fmpr_manref(x));

    if ((rnd == FMPR_RND_DOWN && xsign != sign) ||
        (rnd == FMPR_RND_UP && xsign == sign) ||
        (rnd == FMPR_RND_FLOOR && sign < 0) ||
        (rnd == FMPR_RND_CEIL && sign > 0))
    {
        bits = fmpz_bits(fmpr_manref(x));
        shift = 2 + FLINT_MAX(0, prec - bits);

        fmpz_mul_2exp(fmpr_manref(z), fmpr_manref(x), shift);
        fmpz_sub_ui(fmpr_expref(z), fmpr_expref(x), shift);

        if (sign > 0)
            fmpz_add_ui(fmpr_manref(z), fmpr_manref(z), 1);
        else
            fmpz_sub_ui(fmpr_manref(z), fmpr_manref(z), 1);

        return _fmpr_normalise(fmpr_manref(z), fmpr_expref(z), prec, rnd);
    }
    else
    {
        slong ret = fmpr_set_round(z, x, prec, rnd);

        if (ret == FMPR_RESULT_EXACT)
            return prec - fmpz_bits(fmpr_manref(z));
        else
            return ret;
    }
}
