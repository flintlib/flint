/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "fexpr.h"
#include "fexpr_builtin.h"

void
fexpr_set_d(fexpr_t res, double x)
{
    arf_t t;
    arf_init(t);
    arf_set_d(t, x);  /* no need to free */
    fexpr_set_arf(res, t);
}

void
fexpr_set_re_im_d(fexpr_t res, double x, double y)
{
    if (y == 0.0)
    {
        fexpr_set_d(res, x);
    }
    else if (x == 0.0)
    {
        if (y == 1.0)
        {
            fexpr_set_symbol_builtin(res, FEXPR_NumberI);
        }
        else if (y == -1.0)
        {
            fexpr_set_symbol_builtin(res, FEXPR_NumberI);
            fexpr_neg(res, res);
        }
        else
        {
            fexpr_t im, t;
            fexpr_init(im);
            fexpr_init(t);
            fexpr_set_d(im, y);
            fexpr_set_symbol_builtin(t, FEXPR_NumberI);
            fexpr_mul(res, im, t);
            fexpr_clear(im);
            fexpr_clear(t);
        }
    }
    else
    {
        fexpr_t re, im, t;
        fexpr_init(re);
        fexpr_init(im);
        fexpr_init(t);
        fexpr_set_d(re, x);
        fexpr_set_d(im, y);
        fexpr_set_symbol_builtin(t, FEXPR_NumberI);
        fexpr_mul(res, im, t);
        fexpr_swap(t, res);
        fexpr_add(res, re, t);
        fexpr_clear(re);
        fexpr_clear(im);
        fexpr_clear(t);
    }
}
