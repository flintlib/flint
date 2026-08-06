/*
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"

void
gr_vec_set_shallow(gr_vec_t dest, const gr_vec_t src, gr_ctx_t ctx)
{
    if (dest != src)
    {
        gr_vec_set_length(dest, src->length, ctx);
        _gr_vec_set_shallow(dest->entries, src->entries, src->length, ctx);
    }
}
