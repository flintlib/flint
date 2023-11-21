/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

void
fexpr_call4(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2, const fexpr_t x3, const fexpr_t x4)
{
    slong res_size, f_size, x1_size, x2_size, x3_size, x4_size;
    mp_ptr out;

    f_size = fexpr_size(f);
    x1_size = fexpr_size(x1);
    x2_size = fexpr_size(x2);
    x3_size = fexpr_size(x3);
    x4_size = fexpr_size(x4);

    res_size = FEXPR_HEADER_SIZE + f_size + x1_size + x2_size + x3_size + x4_size;
    fexpr_fit_size(res, res_size);

    out = res->data;
    out[0] = FEXPR_TYPE_CALL4 | (res_size << FEXPR_TYPE_BITS);
    out += FEXPR_HEADER_SIZE;
    flint_mpn_copyi(out, f->data, f_size); out += f_size;
    flint_mpn_copyi(out, x1->data, x1_size); out += x1_size;
    flint_mpn_copyi(out, x2->data, x2_size); out += x2_size;
    flint_mpn_copyi(out, x3->data, x3_size); out += x3_size;
    flint_mpn_copyi(out, x4->data, x4_size);
}
