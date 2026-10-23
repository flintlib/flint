/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fq.h"
#include "fq_poly.h"
#include "gr_poly.h"

void
_fq_poly_mul_classical(fq_struct * rop,
                                  const fq_struct * op1, slong len1,
                                  const fq_struct * op2, slong len2,
                                  const fq_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fq_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_mullow_classical(rop, op1, len1, op2, len2, len1 + len2 - 1, gr_ctx));
}

void
fq_poly_mul_classical(fq_poly_t rop, const fq_poly_t op1,
                                       const fq_poly_t op2, const fq_ctx_t ctx)
{
    const slong len = op1->length + op2->length - 1;

    if (op1->length == 0 || op2->length == 0)
    {
        fq_poly_zero(rop, ctx);
        return;
    }

    if (rop == op1 || rop == op2)
    {
        fq_poly_t t;

        fq_poly_init2(t, len, ctx);
        _fq_poly_mul_classical(t->coeffs, op1->coeffs, op1->length,
                                          op2->coeffs, op2->length, ctx);
        fq_poly_swap(rop, t, ctx);
        fq_poly_clear(t, ctx);
    }
    else
    {
        fq_poly_fit_length(rop, len, ctx);
        _fq_poly_mul_classical(rop->coeffs, op1->coeffs,
                                          op1->length, op2->coeffs,
                                          op2->length, ctx);
    }

    _fq_poly_set_length(rop, len, ctx);
}
