/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

static __inline__ void
fmpz_addmul_2exp(fmpz_t z, const fmpz_t x, const fmpz_t y, ulong shift)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mul_2exp(t, y, shift);
    fmpz_add(z, x, t);
    fmpz_clear(t);
}

static slong _fmpr_add_special(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_zero(x))
    {
        if (fmpr_is_zero(y))
        {
            fmpr_zero(z);
            return FMPR_RESULT_EXACT;
        }
        else
            return fmpr_set_round(z, y, prec, rnd);
    }
    else if (fmpr_is_zero(y))
    {
        return fmpr_set_round(z, x, prec, rnd);
    }
    else if (fmpr_is_nan(x) || fmpr_is_nan(y)
        || (fmpr_is_pos_inf(x) && fmpr_is_neg_inf(y))
        || (fmpr_is_neg_inf(x) && fmpr_is_pos_inf(y)))
    {
        fmpr_nan(z);
        return FMPR_RESULT_EXACT;
    }
    else if (fmpr_is_special(x))
    {
        fmpr_set(z, x);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        fmpr_set(z, y);
        return FMPR_RESULT_EXACT;
    }
}

slong
fmpr_add_naive(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)
{
    slong shift, xsize, ysize;

    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        return _fmpr_add_special(z, x, y, prec, rnd);
    }

    shift = _fmpz_sub_small(fmpr_expref(x), fmpr_expref(y));

    if (shift == 0)
    {
        fmpz_add(fmpr_manref(z), fmpr_manref(x), fmpr_manref(y));
        fmpz_set(fmpr_expref(z), fmpr_expref(x));
    }
    else if (shift > 0)
    {
        ysize = _fmpz_size(fmpr_manref(y)) * FLINT_BITS;

        /* x and y do not overlap */
        if (shift > ysize && prec != FMPR_PREC_EXACT)
        {
            /* y does not overlap with result */
            if (ysize + prec - (slong) fmpz_bits(fmpr_manref(x)) < shift)
            {
                return _fmpr_add_eps(z, x, fmpz_sgn(fmpr_manref(y)), prec, rnd);
            }
        }

        fmpz_addmul_2exp(fmpr_manref(z), fmpr_manref(y), fmpr_manref(x), shift);
        fmpz_set(fmpr_expref(z), fmpr_expref(y));
    }
    else
    {
        shift = -shift;

        xsize = _fmpz_size(fmpr_manref(x)) * FLINT_BITS;

        /* x and y do not overlap */
        if (shift > xsize && prec != FMPR_PREC_EXACT)
        {
            /* y does not overlap with result */
            if (xsize + prec - (slong) fmpz_bits(fmpr_manref(y)) < shift)
            {
                return _fmpr_add_eps(z, y, fmpz_sgn(fmpr_manref(x)), prec, rnd);
            }
        }

        fmpz_addmul_2exp(fmpr_manref(z), fmpr_manref(x), fmpr_manref(y), shift);
        fmpz_set(fmpr_expref(z), fmpr_expref(x));
    }

    return _fmpr_normalise(fmpr_manref(z), fmpr_expref(z), prec, rnd);
}

