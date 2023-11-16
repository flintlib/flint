/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

static int
use_algebraic(const fmpz_t v, const fmpz_t w, slong prec)
{
    fmpz q = *w;
    int r;

    if (COEFF_IS_MPZ(q))
        return 0;

    if (q <= 6)
        return 1;

    r = flint_ctz(q);
    q >>= r;

    if (r >= 4 && prec < (r - 3) * 300)
        return 0;

    if (q > 1000)
        return 0;

    if (prec < 1500 + 150 * q)
        return 0;

    return 1;
}

void
_arb_sin_cos_pi_fmpq_oct(arb_t s, arb_t c,
            const fmpz_t v, const fmpz_t w, slong prec)
{
    if (use_algebraic(v, w, prec))
    {
        _arb_sin_cos_pi_fmpq_algebraic(s, c, *v, *w, prec);
    }
    else
    {
        arb_const_pi(s, prec);
        arb_mul_fmpz(s, s, v, prec);
        arb_div_fmpz(s, s, w, prec);
        arb_sin_cos(s, c, s, prec);
    }
}

void
_arb_sin_pi_fmpq_oct(arb_t s, const fmpz_t v, const fmpz_t w, slong prec)
{
    if (use_algebraic(v, w, prec))
    {
        _arb_sin_pi_fmpq_algebraic(s, *v, *w, prec);
    }
    else
    {
        arb_const_pi(s, prec);
        arb_mul_fmpz(s, s, v, prec);
        arb_div_fmpz(s, s, w, prec);
        arb_sin(s, s, prec);
    }
}

void
_arb_cos_pi_fmpq_oct(arb_t c, const fmpz_t v, const fmpz_t w, slong prec)
{
    if (use_algebraic(v, w, prec))
    {
        _arb_cos_pi_fmpq_algebraic(c, *v, *w, prec);
    }
    else
    {
        arb_const_pi(c, prec);
        arb_mul_fmpz(c, c, v, prec);
        arb_div_fmpz(c, c, w, prec);
        arb_cos(c, c, prec);
    }
}

static unsigned int
reduce_octant(fmpz_t v, fmpz_t w, const fmpq_t x)
{
    const fmpz * p = fmpq_numref(x);
    const fmpz * q = fmpq_denref(x);
    unsigned int octant;
    flint_bitcnt_t vval, wval;

    if (*p > COEFF_MIN / 8 &&
        *p < COEFF_MAX / 8 &&
        *q > 0             &&
        *q < COEFF_MAX / 4)
    {
        slong pp, qq, ww, vv, tt;

        pp = *p;
        qq = *q;

        tt = pp * 4;
        ww = tt / qq;
        vv = tt - qq * ww;
        /* compute correct (floor) quotient and remainder */
        if (vv < 0)
        {
            ww--;
            vv += qq;
        }

        octant = ((ulong) ww) % 8;
        ww = qq * 4;
        if (octant % 2 != 0)
            vv = qq - vv;

        if (vv != 0)
        {
            vval = flint_ctz(vv);
            wval = flint_ctz(ww);
            vval = FLINT_MIN(vval, wval);
            vv >>= vval;
            ww >>= vval;
        }

        fmpz_set_si(v, vv);
        fmpz_set_si(w, ww);
    }
    else
    {
        fmpz_mul_2exp(w, p, 2);
        fmpz_fdiv_qr(w, v, w, q);
        octant = fmpz_fdiv_ui(w, 8);
        fmpz_mul_2exp(w, q, 2);

        if (octant % 2 != 0)
            fmpz_sub(v, q, v);

        vval = fmpz_val2(v);
        wval = fmpz_val2(w);
        vval = FLINT_MIN(vval, wval);

        if (vval != 0)
        {
            fmpz_tdiv_q_2exp(v, v, vval);
            fmpz_tdiv_q_2exp(w, w, vval);
        }
    }

    return octant;
}

void
arb_sin_cos_pi_fmpq(arb_t s, arb_t c, const fmpq_t x, slong prec)
{
    fmpz_t v, w;
    unsigned int octant;

    fmpz_init(v);
    fmpz_init(w);

    octant = reduce_octant(v, w, x);

    if ((octant + 1) % 4 < 2)
        _arb_sin_cos_pi_fmpq_oct(s, c, v, w, prec);
    else
        _arb_sin_cos_pi_fmpq_oct(c, s, v, w, prec);

    if ((octant + 6) % 8 < 4)  arb_neg(c, c);
    if (octant >= 4)           arb_neg(s, s);

    fmpz_clear(v);
    fmpz_clear(w);
}

void
arb_sin_pi_fmpq(arb_t s, const fmpq_t x, slong prec)
{
    fmpz_t v, w;
    unsigned int octant;

    fmpz_init(v);
    fmpz_init(w);

    octant = reduce_octant(v, w, x);

    if ((octant + 1) % 4 < 2)
        _arb_sin_pi_fmpq_oct(s, v, w, prec);
    else
        _arb_cos_pi_fmpq_oct(s, v, w, prec);

    if (octant >= 4)
        arb_neg(s, s);

    fmpz_clear(v);
    fmpz_clear(w);
}

void
arb_cos_pi_fmpq(arb_t c, const fmpq_t x, slong prec)
{
    fmpz_t v, w;
    unsigned int octant;

    fmpz_init(v);
    fmpz_init(w);

    octant = reduce_octant(v, w, x);

    if ((octant + 1) % 4 < 2)
        _arb_cos_pi_fmpq_oct(c, v, w, prec);
    else
        _arb_sin_pi_fmpq_oct(c, v, w, prec);

    if ((octant + 6) % 8 < 4)
        arb_neg(c, c);

    fmpz_clear(v);
    fmpz_clear(w);
}

