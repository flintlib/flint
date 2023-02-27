/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

/* TODO: correct rounding with negative b */

void
fmpr_pow_sloppy_fmpz(fmpr_t y, const fmpr_t b, const fmpz_t e,
    slong prec, fmpr_rnd_t rnd)
{
    slong i, wp, bits;

    if (fmpz_is_zero(e))
    {
        fmpr_set_ui(y, UWORD(1));
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);
        fmpr_pow_sloppy_fmpz(y, b, f, prec + 2,
            (rnd == FMPR_RND_FLOOR || rnd == FMPR_RND_DOWN)
            ? FMPR_RND_UP : FMPR_RND_DOWN);
        fmpr_ui_div(y, UWORD(1), y, prec, rnd);
        fmpz_clear(f);
        return;
    }

    if (y == b)
    {
        fmpr_t t;
        fmpr_init(t);
        fmpr_set(t, b);
        fmpr_pow_sloppy_fmpz(y, t, e, prec, rnd);
        fmpr_clear(t);
        return;
    }

    fmpr_set(y, b);

    bits = fmpz_bits(e);
    wp = FMPR_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        fmpr_mul(y, y, y, wp, rnd);
        if (fmpz_tstbit(e, i))
            fmpr_mul(y, y, b, wp, rnd);
    }
}

void
fmpr_pow_sloppy_ui(fmpr_t y, const fmpr_t b, ulong e, slong prec, fmpr_rnd_t rnd)
{
    fmpz_t f;
    fmpz_init_set_ui(f, e);
    fmpr_pow_sloppy_fmpz(y, b, f, prec, rnd);
    fmpz_clear(f);
}

void
fmpr_pow_sloppy_si(fmpr_t y, const fmpr_t b, slong e, slong prec, fmpr_rnd_t rnd)
{
    fmpz_t f;
    fmpz_init(f);
    fmpz_set_si(f, e);
    fmpr_pow_sloppy_fmpz(y, b, f, prec, rnd);
    fmpz_clear(f);
}
