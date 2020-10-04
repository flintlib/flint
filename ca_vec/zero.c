/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_vec.h"

void
_ca_vec_zero(ca_ptr res, slong len, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        ca_zero(res + i, ctx);
}

void
ca_vec_zero(ca_vec_t res, slong len, ca_ctx_t ctx)
{
    ca_vec_set_length(res, len, ctx);
    _ca_vec_zero(ca_vec_entry(res, 0), ca_vec_length(res, ctx), ctx);
}
