/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

static void
_fmpr_div_special(fmpr_t z, const fmpr_t x, const fmpr_t y)
{
    if ((fmpr_is_zero(x) && !fmpr_is_zero(y) && !fmpr_is_nan(y)) ||
        (fmpr_is_inf(y) && !fmpr_is_special(x)))
    {
        fmpr_zero(z);
    }
    else if (fmpr_is_zero(y) || (fmpr_is_special(x) && fmpr_is_special(y)) ||
        fmpr_is_nan(x) || fmpr_is_nan(y))
    {
        fmpr_nan(z);
    }
    else if (fmpr_sgn(x) == fmpr_sgn(y))
        fmpr_pos_inf(z);
    else
        fmpr_neg_inf(z);
}

static __inline__ int
low_bits_are_zero(const fmpz_t u, int bits)
{
    fmpz f = *u;
    mp_limb_t low;

    if (!COEFF_IS_MPZ(f))
        low = FLINT_ABS(f);
    else
        low = COEFF_TO_PTR(f)->_mp_d[0];

    return (low & ((UWORD(1) << bits) - 1)) == 0;
}

slong
fmpr_div(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        _fmpr_div_special(z, x, y);
        return FMPR_RESULT_EXACT;
    }

    /* division by power of two <=> shift exponents */
    if (fmpz_is_pm1(fmpr_manref(y)))
    {
        if (fmpz_is_one(fmpr_manref(y)))
            fmpz_set(fmpr_manref(z), fmpr_manref(x));
        else
            fmpz_neg(fmpr_manref(z), fmpr_manref(x));
        fmpz_sub(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y));
        return _fmpr_normalise(fmpr_manref(z), fmpr_expref(z), prec, rnd);
    }
    else
    {
        slong xbits, ybits, extra, extra_pad, extra_control;
        int negative;
        fmpz_t t, u;

        /* todo: work out exact needed shift */
        xbits = fmpz_bits(fmpr_manref(x));
        ybits = fmpz_bits(fmpr_manref(y));

        extra = prec - xbits + ybits;
        extra = FLINT_MAX(extra, 0);

        extra_pad = 32;
        extra_control = 24;
        extra += extra_pad;

        fmpz_init(t);
        fmpz_init(u);

        fmpz_mul_2exp(t, fmpr_manref(x), extra);
        fmpz_tdiv_q(u, t, fmpr_manref(y));

        if (low_bits_are_zero(u, extra_control))
        {
            fmpz_t v;
            fmpz_init(v);
            fmpz_mul(v, u, fmpr_manref(y));

            negative = fmpz_sgn(fmpr_manref(x)) != fmpz_sgn(fmpr_manref(y));

            if (!fmpz_equal(t, v))
            {
                if (negative)
                    fmpz_sub_ui(u, u, 1);
                else
                    fmpz_add_ui(u, u, 1);
            }

            fmpz_clear(v);
        }

        fmpz_swap(fmpr_manref(z), u);

        fmpz_clear(t);
        fmpz_clear(u);

        fmpz_sub(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y));
        fmpz_sub_ui(fmpr_expref(z), fmpr_expref(z), extra);

        return _fmpr_normalise(fmpr_manref(z), fmpr_expref(z), prec, rnd);
    }
}

slong
fmpr_div_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd)
{
    fmpr_t t; slong r;
    fmpr_init(t);
    fmpr_set_ui(t, y);
    r = fmpr_div(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

slong
fmpr_ui_div(fmpr_t z, ulong x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)
{
    fmpr_t t; slong r;
    fmpr_init(t);
    fmpr_set_ui(t, x);
    r = fmpr_div(z, t, y, prec, rnd);
    fmpr_clear(t);
    return r;
}

slong
fmpr_div_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd)
{
    fmpr_t t; slong r;
    fmpr_init(t);
    fmpr_set_si(t, y);
    r = fmpr_div(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

slong
fmpr_si_div(fmpr_t z, slong x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)
{
    fmpr_t t; slong r;
    fmpr_init(t);
    fmpr_set_si(t, x);
    r = fmpr_div(z, t, y, prec, rnd);
    fmpr_clear(t);
    return r;
}

slong
fmpr_div_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd)
{
    fmpr_t t; slong r;
    fmpr_init(t);
    fmpr_set_fmpz(t, y);
    r = fmpr_div(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

slong
fmpr_fmpz_div(fmpr_t z, const fmpz_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)
{
    fmpr_t t; slong r;
    fmpr_init(t);
    fmpr_set_fmpz(t, x);
    r = fmpr_div(z, t, y, prec, rnd);
    fmpr_clear(t);
    return r;
}

slong
fmpr_fmpz_div_fmpz(fmpr_t z, const fmpz_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd)
{
    fmpr_t t, u; slong r;
    fmpr_init(t);
    fmpr_init(u);
    fmpr_set_fmpz(t, x);
    fmpr_set_fmpz(u, y);
    r = fmpr_div(z, t, u, prec, rnd);
    fmpr_clear(t);
    fmpr_clear(u);
    return r;
}
