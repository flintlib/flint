/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"

void gr_vec_clear(gr_vec_t vec, gr_ctx_t ctx)
{
    if (vec->alloc != 0)
    {
        _gr_vec_clear(vec->entries, vec->alloc, ctx);
        flint_free(vec->entries);
    }
}
