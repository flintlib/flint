/*
   Copyright (C) 2025 Marc Mezzarobba

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
   */

#include "gr_vec.h"
#include "perm.h"

void
_gr_vec_permute(gr_ptr vec, slong * perm, slong len, gr_ctx_t ctx)
{
    for (slong base = 0; base < len; base++)
    {
        slong k = base;
        while (perm[k] != k)
        {
            gr_swap(GR_ENTRY(vec, base, ctx->sizeof_elem),
                    GR_ENTRY(vec, k, ctx->sizeof_elem),
                    ctx);
            FLINT_SWAP(slong, perm[k], k);
        }
    }
}

int
gr_vec_permute(gr_vec_t dest, gr_vec_t src, slong * perm, gr_ctx_t ctx)
{
    if (dest != src)
    {
        int status = gr_vec_set(dest, src, ctx);
        if (status != GR_SUCCESS)
            return status;
    }

    slong * _perm = _perm_init(src->length);
    _perm_set(_perm, perm, src->length);

    _gr_vec_permute(dest->entries, _perm, dest->length, ctx);

    _perm_clear(_perm);

    return GR_SUCCESS;
}
