/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_sqr_classical) (TEMPLATE(T, struct) * rop,
                                  const TEMPLATE(T, struct) * op, slong len,
                                  const TEMPLATE(T, ctx_t) ctx)
{
    if (len == 1)
    {
        TEMPLATE(T, mul) (rop, op, op, ctx);
    }
    else
    {
        slong i;
        TEMPLATE(T, t) t;

        TEMPLATE(T, init) (t, ctx);

        _TEMPLATE(T, TEMPLATE(poly_scalar_mul, T)) (rop, op, len, op, ctx);

        _TEMPLATE(T, TEMPLATE(poly_scalar_mul, T)) (rop + len, op + 1, len - 1,
                                                    op + len - 1, ctx);

        for (i = 1; i < len - 1; i++)
            _TEMPLATE(T, TEMPLATE(poly_scalar_addmul, T)) (rop + i + 1, op + 1,
                                                           i - 1, op + i, ctx);

        for (i = 1; i < 2 * len - 2; i++)
            TEMPLATE(T, add) (rop + i, rop + i, rop + i, ctx);

        for (i = 1; i < len - 1; i++)
        {
            TEMPLATE(T, sqr) (t, op + i, ctx);
            TEMPLATE(T, add) (rop + 2 * i, rop + 2 * i, t, ctx);
        }
        TEMPLATE(T, clear) (t, ctx);
    }
}

void
TEMPLATE(T, poly_sqr_classical) (TEMPLATE(T, poly_t) rop,
                                 const TEMPLATE(T, poly_t) op,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    const slong len = 2 * op->length - 1;

    if (op->length == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
        return;
    }

    if (rop == op)
    {
        TEMPLATE(T, poly_t) t;

        TEMPLATE(T, poly_init2) (t, len, ctx);
        _TEMPLATE(T, poly_sqr_classical) (t->coeffs, op->coeffs, op->length,
                                          ctx);
        TEMPLATE(T, poly_swap) (rop, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, len, ctx);
        _TEMPLATE(T, poly_sqr_classical) (rop->coeffs, op->coeffs, op->length,
                                          ctx);
    }

    _TEMPLATE(T, poly_set_length) (rop, len, ctx);
}


#endif
