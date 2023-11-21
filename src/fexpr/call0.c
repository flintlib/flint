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
fexpr_call0(fexpr_t res, const fexpr_t f)
{
    slong res_size, f_size;
    mp_ptr out;

    f_size = fexpr_size(f);

    res_size = FEXPR_HEADER_SIZE + f_size;
    fexpr_fit_size(res, res_size);

    out = res->data;
    out[0] = FEXPR_TYPE_CALL0 | (res_size << FEXPR_TYPE_BITS);
    out += FEXPR_HEADER_SIZE;
    flint_mpn_copyi(out, f->data, f_size);
}
