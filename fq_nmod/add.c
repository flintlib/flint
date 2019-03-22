/*
    Copyright (C) 2011, 2012 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void fq_nmod_add(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
{
    slong max = FLINT_MAX(op1->length, op2->length);

    nmod_polydr_fit_length(rop, max, ctx->fpctx);

    _nmod_poly_add(rop->coeffs, 
                   op1->coeffs, op1->length, op2->coeffs, op2->length, 
                   ctx->fpctx->mod);

    _nmod_polydr_set_length(rop, max);
    _nmod_polydr_normalise(rop);
}
