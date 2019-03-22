/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void fq_nmod_neg(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    slong len = op->length;

    nmod_polydr_fit_length(rop, len, ctx->fpctx);
    _nmod_polydr_set_length(rop, len);

    nmod_polydr_neg(rop, op, ctx->fpctx);
}
