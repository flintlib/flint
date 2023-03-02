/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "fmpz_poly.h"

void
_TEMPLATE(T, poly_sqr_KS) (TEMPLATE(T, struct) * rop,
                           const TEMPLATE(T, struct) * op, slong len,
                           const TEMPLATE(T, ctx_t) ctx)
{
    const slong in_len = len;
    const slong d = TEMPLATE(T, ctx_degree) (ctx);
    slong bits, i;
    fmpz *f, *g;

    TEMPLATE(CAP_T, VEC_NORM) (op, len, ctx);

    if (!len)
    {
        if (2 * in_len - 1 > 0)
            _TEMPLATE(T, poly_zero) (rop, 2 * in_len - 1, ctx);
        return;
    }

    bits = 2 * fmpz_bits(TEMPLATE(T, ctx_prime) (ctx))
        + FLINT_BIT_COUNT(d) + FLINT_BIT_COUNT(len);

    f = _fmpz_vec_init((2 * len - 1) + len);
    g = f + (2 * len - 1);

    for (i = 0; i < len; i++)
    {
        TEMPLATE(T, bit_pack) (g + i, op + i, bits, ctx);
    }

    _fmpz_poly_sqr(f, g, len);

    for (i = 0; i < 2 * len - 1; i++)
    {
        TEMPLATE(T, bit_unpack) (rop + i, f + i, bits, ctx);
    }

    _TEMPLATE(T, poly_zero) (rop + (2 * len - 1), 2 * (in_len - len), ctx);

    _fmpz_vec_clear(f, (2 * len - 1) + len);
}

void
TEMPLATE(T, poly_sqr_KS) (TEMPLATE(T, poly_t) rop,
                          const TEMPLATE(T, poly_t) op,
                          const TEMPLATE(T, ctx_t) ctx)
{
    const slong len = 2 * op->length - 1;

    if (op->length == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, len, ctx);
        _TEMPLATE(T, poly_sqr_KS) (rop->coeffs, op->coeffs, op->length, ctx);
        _TEMPLATE(T, poly_set_length) (rop, len, ctx);
    }
}


#endif
