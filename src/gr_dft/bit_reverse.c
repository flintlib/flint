/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr.h"
#include "gr_dft.h"

void
_gr_dft_bit_reverse(gr_ptr x, slong stride, int depth, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    ulong i, j, n = UWORD(1) << depth;

    for (i = 0; i < n; i++)
    {
        j = n_revbin(i, depth);

        if (i < j)
            gr_swap(GR_ENTRY(x, (slong) i * stride, sz),
                    GR_ENTRY(x, (slong) j * stride, sz), ctx);
    }
}
