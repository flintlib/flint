/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acf.h"

void
acf_approx_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_zero(acf_imagref(x)))
    {
        /* todo: arf_inv */
        arf_ui_div(acf_realref(res), 1, acf_realref(x), prec, rnd);
        arf_zero(acf_imagref(res));
    }
    else if (arf_is_zero(acf_realref(x)))
    {
        /* todo: arf_inv */
        arf_si_div(acf_imagref(res), -1, acf_imagref(x), prec, rnd);
        arf_zero(acf_realref(res));
    }
    else
    {
        arf_t t;
        arf_init(t);

        arf_sosq(t, acf_realref(x), acf_imagref(x), prec, rnd);
        arf_div(acf_realref(res), acf_realref(x), t, prec, rnd);
        arf_div(acf_imagref(res), acf_imagref(x), t, prec, rnd);
        arf_neg(acf_imagref(res), acf_imagref(res));

        arf_clear(t);
    }
}

#define a acf_realref(x)
#define b acf_imagref(x)
#define c acf_realref(y)
#define d acf_imagref(y)

void
acf_approx_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    if (arf_is_zero(d))
    {
        if (arf_is_zero(b))
        {
            arf_div(acf_realref(res), a, c, prec, rnd);
            arf_zero(acf_imagref(res));
        }
        else if (arf_is_zero(a))
        {
            arf_div(acf_imagref(res), b, c, prec, rnd);
            arf_zero(acf_realref(res));
        }
        else if (res != y)
        {
            arf_div(acf_realref(res), a, c, prec, rnd);
            arf_div(acf_imagref(res), b, c, prec, rnd);
        }
        else
        {
            arf_t t;
            arf_init(t);
            arf_set(t, c);
            arf_div(acf_realref(res), a, t, prec, rnd);
            arf_div(acf_imagref(res), b, t, prec, rnd);
            arf_clear(t);
        }
    }
    else if (arf_is_zero(c))
    {
        if (arf_is_zero(b))
        {
            arf_div(acf_imagref(res), a, d, prec, rnd);
            arf_neg(acf_imagref(res), acf_imagref(res));
            arf_zero(acf_realref(res));
        }
        else if (arf_is_zero(a))
        {
            arf_div(acf_realref(res), b, d, prec, rnd);
            arf_zero(acf_imagref(res));
        }
        else if (res != y)
        {
            arf_div(acf_realref(res), a, d, prec, rnd);
            arf_div(acf_imagref(res), b, d, prec, rnd);
            arf_swap(acf_realref(res), acf_imagref(res));
            arf_neg(acf_imagref(res), acf_imagref(res));
        }
        else
        {
            arf_t t;
            arf_init(t);
            arf_set(t, d);
            arf_div(acf_realref(res), a, t, prec, rnd);
            arf_div(acf_imagref(res), b, t, prec, rnd);
            arf_swap(acf_realref(res), acf_imagref(res));
            arf_neg(acf_imagref(res), acf_imagref(res));
            arf_clear(t);
        }
    }
    else
    {
        arf_t t;
        acf_t u;

        arf_init(t);

        arf_sosq(t, acf_realref(y), acf_imagref(y), prec, rnd);

        arf_init_set_shallow(acf_realref(u), acf_realref(y));
        arf_init_neg_shallow(acf_imagref(u), acf_imagref(y));

        acf_mul(res, x, u, prec, rnd);

        arf_div(acf_realref(res), acf_realref(res), t, prec, rnd);
        arf_div(acf_imagref(res), acf_imagref(res), t, prec, rnd);

        arf_clear(t);
    }
}

#undef a
#undef b
#undef c
#undef d
