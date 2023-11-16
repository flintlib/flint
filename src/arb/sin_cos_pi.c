/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_sin_pi(arb_t y, const arb_t x, slong prec)
{
    arb_t t;
    arb_t u;
    fmpz_t v;

    if (!arb_is_finite(x))
    {
        arb_indeterminate(y);
        return;
    }

    if (arf_cmpabs_2exp_si(arb_midref(x), FLINT_MAX(65536, (4*prec))) > 0)
    {
        arf_zero(arb_midref(y));
        mag_one(arb_radref(y));
        return;
    }

    arb_init(t);
    arb_init(u);
    fmpz_init(v);

    arb_mul_2exp_si(t, x, 1);
    arf_get_fmpz(v, arb_midref(t), ARF_RND_NEAR);
    arb_sub_fmpz(t, t, v, prec);

    arb_const_pi(u, prec);
    arb_mul(t, t, u, prec);
    arb_mul_2exp_si(t, t, -1);

    switch (fmpz_fdiv_ui(v, 4))
    {
        case 0:
            arb_sin(y, t, prec);
            break;
        case 1:
            arb_cos(y, t, prec);
            break;
        case 2:
            arb_sin(y, t, prec);
            arb_neg(y, y);
            break;
        default:
            arb_cos(y, t, prec);
            arb_neg(y, y);
            break;
    }

    fmpz_clear(v);
    arb_clear(t);
    arb_clear(u);
}

void
arb_cos_pi(arb_t y, const arb_t x, slong prec)
{
    arb_t t;
    arb_t u;
    fmpz_t v;

    if (!arb_is_finite(x))
    {
        arb_indeterminate(y);
        return;
    }

    if (arf_cmpabs_2exp_si(arb_midref(x), FLINT_MAX(65536, (4*prec))) > 0)
    {
        arf_zero(arb_midref(y));
        mag_one(arb_radref(y));
        return;
    }

    arb_init(t);
    arb_init(u);
    fmpz_init(v);

    arb_mul_2exp_si(t, x, 1);
    arf_get_fmpz(v, arb_midref(t), ARF_RND_NEAR);
    arb_sub_fmpz(t, t, v, prec);

    arb_const_pi(u, prec);
    arb_mul(t, t, u, prec);
    arb_mul_2exp_si(t, t, -1);

    switch (fmpz_fdiv_ui(v, 4))
    {
        case 0:
            arb_cos(y, t, prec);
            break;
        case 1:
            arb_sin(y, t, prec);
            arb_neg(y, y);
            break;
        case 2:
            arb_cos(y, t, prec);
            arb_neg(y, y);
            break;
        default:
            arb_sin(y, t, prec);
            break;
    }

    fmpz_clear(v);
    arb_clear(t);
    arb_clear(u);
}

void
arb_sin_cos_pi(arb_t s, arb_t c, const arb_t x, slong prec)
{
    arb_t t;
    arb_t u;
    fmpz_t v;

    if (!arb_is_finite(x))
    {
        arb_indeterminate(s);
        arb_indeterminate(c);
        return;
    }

    if (arf_cmpabs_2exp_si(arb_midref(x), FLINT_MAX(65536, (4*prec))) > 0)
    {
        arf_zero(arb_midref(s));
        mag_one(arb_radref(s));
        arf_zero(arb_midref(c));
        mag_one(arb_radref(c));
        return;
    }

    arb_init(t);
    arb_init(u);
    fmpz_init(v);

    arb_mul_2exp_si(t, x, 1);
    arf_get_fmpz(v, arb_midref(t), ARF_RND_NEAR);
    arb_sub_fmpz(t, t, v, prec);

    arb_const_pi(u, prec);
    arb_mul(t, t, u, prec);
    arb_mul_2exp_si(t, t, -1);

    switch (fmpz_fdiv_ui(v, 4))
    {
        case 0:
            arb_sin_cos(s, c, t, prec);
            break;
        case 1:
            arb_sin_cos(c, s, t, prec);
            arb_neg(c, c);
            break;
        case 2:
            arb_sin_cos(s, c, t, prec);
            arb_neg(s, s);
            arb_neg(c, c);
            break;
        default:
            arb_sin_cos(c, s, t, prec);
            arb_neg(s, s);
            break;
    }

    fmpz_clear(v);
    arb_clear(t);
    arb_clear(u);
}

