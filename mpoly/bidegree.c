/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* deg_gen(0) and deg_gen(1) of lt(A) packed into one ulong */
ulong _mpoly_bidegree(
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx)
{
    slong off0, shift0, off1, shift1;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);

    FLINT_ASSERT(mctx->ord == ORD_LEX);

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, Abits, mctx);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, Abits, mctx);

    return pack_exp2((Aexps[off0] >> shift0) & mask,
                     (Aexps[off1] >> shift1) & mask);
}
