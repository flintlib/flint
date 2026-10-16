/*
    Copyright (C) 2010, 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "gr_poly.h"

void
_TEMPLATE(T, poly_pow) (TEMPLATE(T, struct) * rop,
                        const TEMPLATE(T, struct) * op, slong len, ulong e,
                        const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_pow_ui_binexp(rop, op, len, e, gr_ctx));
}

void
TEMPLATE(T, poly_pow) (TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                       ulong e, const TEMPLATE(T, ctx_t) ctx)
{
    const slong len = op->length;

    if ((len < 2) | (e < UWORD(3)))
    {
        if (e == UWORD(0))
            TEMPLATE(T, poly_one) (rop, ctx);
        else if (len == 0)
            TEMPLATE(T, poly_zero) (rop, ctx);
        else if (len == 1)
        {
            fmpz_t f;
            fmpz_init_set_ui(f, e);

            TEMPLATE(T, poly_fit_length) (rop, 1, ctx);
            TEMPLATE(T, pow) (rop->coeffs + 0, op->coeffs + 0, f, ctx);
            _TEMPLATE(T, poly_set_length) (rop, 1, ctx);

            fmpz_clear(f);
        }
        else if (e == UWORD(1))
            TEMPLATE(T, poly_set) (rop, op, ctx);
        else                    /* e == UWORD(2) */
            TEMPLATE(T, poly_sqr) (rop, op, ctx);
    }
    else
    {
        const slong rlen = (slong) e * (len - 1) + 1;

        if (rop != op)
        {
            TEMPLATE(T, poly_fit_length) (rop, rlen, ctx);
            _TEMPLATE(T, poly_pow) (rop->coeffs, op->coeffs, len, e, ctx);
            _TEMPLATE(T, poly_set_length) (rop, rlen, ctx);
        }
        else
        {
            TEMPLATE(T, poly_t) t;
            TEMPLATE(T, poly_init2) (t, rlen, ctx);
            _TEMPLATE(T, poly_pow) (t->coeffs, op->coeffs, len, e, ctx);
            _TEMPLATE(T, poly_set_length) (t, rlen, ctx);
            TEMPLATE(T, poly_swap) (rop, t, ctx);
            TEMPLATE(T, poly_clear) (t, ctx);
        }
    }
}


#endif
