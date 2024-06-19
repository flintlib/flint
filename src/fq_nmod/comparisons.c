/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fq_nmod.h"

int fq_nmod_equal(const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t FLINT_UNUSED(ctx))
{
    return nmod_poly_equal(op1, op2);
}

int fq_nmod_is_zero(const fq_nmod_t op, const fq_nmod_ctx_t FLINT_UNUSED(ctx))
{
    return nmod_poly_is_zero(op);
}

int fq_nmod_is_one(const fq_nmod_t op, const fq_nmod_ctx_t FLINT_UNUSED(ctx))
{
    return nmod_poly_is_one(op);
}

int fq_nmod_cmp(const fq_nmod_t a, const fq_nmod_t b, const fq_nmod_ctx_t FLINT_UNUSED(ctx))
{
    slong i;

    if (a->length != b->length)
        return a->length < b->length ? -1 : 1;

    for (i = 0; i < a->length; i++)
    {
        if (a->coeffs[i] != b->coeffs[i])
            return a->coeffs[i] < b->coeffs[i] ? -1 : 1;
    }

    return 0;
}
