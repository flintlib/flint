/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

slong
fmpr_root(fmpr_t y, const fmpr_t x, ulong k, slong prec, fmpr_rnd_t rnd)
{
    slong r;

    if (k == 0)
    {
        fmpr_nan(y);
        return FMPR_RESULT_EXACT;
    }
    else if (k == 1)
    {
        return fmpr_set_round(y, x, prec, rnd);
    }
    else if (k == 2)
    {
        return fmpr_sqrt(y, x, prec, rnd);
    }

    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_zero(y);
        else if (fmpr_is_pos_inf(x))
            fmpr_pos_inf(y);
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }

    if (fmpr_sgn(x) < 0)
    {
        fmpr_nan(y);
        return FMPR_RESULT_EXACT;
    }

    {
        fmpr_t t;
        fmpz_t a, b;

        fmpr_init(t);
        fmpz_init(a);
        fmpz_init(b);

        /* (m * 2^(aq+b))^(1/q) = (m*2^b)^(1/q) * 2^a */
        fmpz_set_ui(a, k);
        fmpz_fdiv_qr(a, b, fmpr_expref(x), a);

        fmpz_set(fmpr_manref(t), fmpr_manref(x));
        fmpz_set(fmpr_expref(t), b);

#if MPFR_VERSION_MAJOR >= 4
        CALL_MPFR_FUNC_K(r, mpfr_rootn_ui, y, t, k, prec, rnd);
#else
        CALL_MPFR_FUNC_K(r, mpfr_root, y, t, k, prec, rnd);
#endif

        fmpr_mul_2exp_fmpz(y, y, a);

        fmpr_clear(t);
        fmpz_clear(a);
        fmpz_clear(b);

        return r;
    }
}

