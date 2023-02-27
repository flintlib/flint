/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

static __inline__ int sgn(int x)
{
    return (x > 0) - (x < 0);
}

int
fmpr_cmp_2exp_si(const fmpr_t x, slong e)
{
    slong bc;
    int ret;
    fmpz_t t;

    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x)) return -1;
        if (fmpr_is_pos_inf(x)) return 1;
        if (fmpr_is_neg_inf(x)) return -1;
        if (fmpr_is_nan(x)) return 0;
        return -1;
    }

    if (fmpz_is_one(fmpr_manref(x)))
        return sgn(fmpz_cmp_si(fmpr_expref(x), e));

    if (fmpz_sgn(fmpr_manref(x)) < 0)
        return -1;

    bc = fmpz_bits(fmpr_manref(x));

    fmpz_init(t);

    fmpz_add_si_inline(t, fmpr_expref(x), bc);
    fmpz_sub_si_inline(t, t, e);

    ret = (fmpz_sgn(t) <= 0) ? -1 : 1;

    fmpz_clear(t);
    return ret;
}

