/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

void
fq_zech_bit_unpack(fq_zech_t rop, const fmpz_t f, flint_bitcnt_t bit_size,
                   const fq_zech_ctx_t ctx)
{
    fq_nmod_t ropn;
    fq_nmod_init(ropn, ctx->fq_nmod_ctx);
    
    fq_nmod_bit_unpack(ropn, f, bit_size, ctx->fq_nmod_ctx);
    fq_zech_set_fq_nmod(rop, ropn, ctx);

    fq_nmod_clear(ropn, ctx->fq_nmod_ctx);
}
