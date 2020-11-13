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
_ca_vec_swap(ca_ptr vec1, ca_ptr vec2, slong len, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        ca_swap(vec1 + i, vec2 + i, ctx);
}
