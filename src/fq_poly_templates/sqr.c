/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#if defined(FQ_ZECH_POLY_H)
/* Kept out of line so that the fast path for small squares
   in _fq_zech_poly_sqr does not need a stack frame. */
FLINT_STATIC_NOINLINE void
_fq_zech_poly_sqr_large(fq_zech_struct * rop,
                        const fq_zech_struct * op, slong len,
                        const fq_zech_ctx_t ctx)
{
    if (_fq_zech_poly_sqr_want_univariate(len, ctx))
        _fq_zech_poly_mul_univariate(rop, op, len, op, len, ctx);
    else
        _fq_zech_poly_sqr_classical(rop, op, len, ctx);
}
#endif

void
_TEMPLATE(T, poly_sqr) (TEMPLATE(T, struct) * rop,
                        const TEMPLATE(T, struct) * op, slong len,
                        const TEMPLATE(T, ctx_t) ctx)
{
#if defined(FQ_ZECH_POLY_H)
    if (len < 2 * FQ_ZECH_POLY_MUL_UNIVARIATE_MIN_LEN(ctx))
        _fq_zech_poly_sqr_classical(rop, op, len, ctx);
    else
        _fq_zech_poly_sqr_large(rop, op, len, ctx);
#else
    if (len < TEMPLATE(CAP_T, SQR_CLASSICAL_CUTOFF))
    {
        _TEMPLATE(T, poly_sqr_classical) (rop, op, len, ctx);
    }
#ifdef USE_SQR_REORDER
    else if (TEMPLATE(T, ctx_degree) (ctx) < 4)
    {
        _TEMPLATE(T, poly_sqr_reorder) (rop, op, len, ctx);
    }
#endif
    else
    {
        _TEMPLATE(T, poly_sqr_KS) (rop, op, len, ctx);
    }
#endif
}

void
TEMPLATE(T, poly_sqr) (TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                       const TEMPLATE(T, ctx_t) ctx)
{
    const slong rlen = 2 * op->length - 1;

    if (op->length == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
        return;
    }

    if (rop == op)
    {
        TEMPLATE(T, poly_t) t;

        TEMPLATE(T, poly_init2) (t, rlen, ctx);
        _TEMPLATE(T, poly_sqr) (t->coeffs, op->coeffs, op->length, ctx);
        TEMPLATE(T, poly_swap) (rop, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, rlen, ctx);
        _TEMPLATE(T, poly_sqr) (rop->coeffs, op->coeffs, op->length, ctx);
    }

    _TEMPLATE(T, poly_set_length) (rop, rlen, ctx);
}


#endif
