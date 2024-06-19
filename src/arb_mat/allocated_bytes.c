/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_mat.h"

slong arb_mat_allocated_bytes(const arb_mat_t x)
{
    return _arb_vec_allocated_bytes(x->entries, x->r * x->c) + x->r * sizeof(arb_ptr);
}
