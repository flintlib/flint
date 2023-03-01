/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

void
fexpr_call2(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2)
{
    slong res_size, f_size, x1_size, x2_size;
    mp_ptr out;

    f_size = fexpr_size(f);
    x1_size = fexpr_size(x1);
    x2_size = fexpr_size(x2);

    res_size = FEXPR_HEADER_SIZE + f_size + x1_size + x2_size;
    fexpr_fit_size(res, res_size);

    out = res->data;
    out[0] = FEXPR_TYPE_CALL2 | (res_size << FEXPR_TYPE_BITS);
    out += FEXPR_HEADER_SIZE;
    flint_mpn_copyi(out, f->data, f_size); out += f_size;
    flint_mpn_copyi(out, x1->data, x1_size); out += x1_size;
    flint_mpn_copyi(out, x2->data, x2_size);
}

void
fexpr_call_builtin2(fexpr_t res, slong f, const fexpr_t x, const fexpr_t y)
{
    fexpr_t t;
    ulong d;
    t->data = &d;
    t->alloc = 1;
    fexpr_set_symbol_builtin(t, f);

    if (res == x || res == y)
    {
        fexpr_t u;
        fexpr_init(u);
        fexpr_call2(u, t, x, y);
        fexpr_swap(res, u);
        fexpr_clear(u);
    }
    else
    {
        fexpr_call2(res, t, x, y);
    }
}
