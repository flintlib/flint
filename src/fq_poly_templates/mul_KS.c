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
_TEMPLATE(T, poly_mul_KS) (TEMPLATE(T, struct) * rop,
                           const TEMPLATE(T, struct) * op1, slong len1,
                           const TEMPLATE(T, struct) * op2, slong len2,
                           const TEMPLATE(T, ctx_t) ctx)
{
    const slong in1_len = len1, in2_len = len2;
    const slong d = TEMPLATE(T, ctx_degree) (ctx);
    slong bits, i;
    fmpz *f, *g, *h;

    TEMPLATE(CAP_T, VEC_NORM) (op1, len1, ctx);
    TEMPLATE(CAP_T, VEC_NORM) (op2, len2, ctx);

    if (!len1 | !len2)
    {
        if (in1_len + in2_len - 1 > 0)
            _TEMPLATE(T, poly_zero) (rop, in1_len + in2_len - 1, ctx);
        return;
    }

    bits = 2 * fmpz_bits(TEMPLATE(T, ctx_prime) (ctx))
        + FLINT_BIT_COUNT(d) + FLINT_BIT_COUNT(FLINT_MIN(len1, len2));

    f = _fmpz_vec_init((len1 + len2 - 1) + (len1) + (len2));
    g = f + (len1 + len2 - 1);
    h = g + len1;

    for (i = 0; i < len1; i++)
    {
        TEMPLATE(T, bit_pack) (g + i, op1 + i, bits, ctx);
    }
    for (i = 0; i < len2; i++)
    {
        TEMPLATE(T, bit_pack) (h + i, op2 + i, bits, ctx);
    }

    if (len1 >= len2)
        _fmpz_poly_mul(f, g, len1, h, len2);
    else
        _fmpz_poly_mul(f, h, len2, g, len1);

    for (i = 0; i < len1 + len2 - 1; i++)
    {
        TEMPLATE(T, bit_unpack) (rop + i, f + i, bits, ctx);
    }

    _TEMPLATE(T, poly_zero) (rop + (len1 + len2 - 1),
                             (in1_len - len1) + (in2_len - len2), ctx);

    _fmpz_vec_clear(f, (len1 + len2 - 1) + (len1) + (len2));
}

void
TEMPLATE(T, poly_mul_KS) (TEMPLATE(T, poly_t) rop,
                          const TEMPLATE(T, poly_t) op1,
                          const TEMPLATE(T, poly_t) op2,
                          const TEMPLATE(T, ctx_t) ctx)
{
    const slong len1 = op1->length;
    const slong len2 = op2->length;
    const slong rlen = len1 + len2 - 1;

    if (len1 == 0 || len2 == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, rlen, ctx);
        _TEMPLATE(T, poly_mul_KS) (rop->coeffs, op1->coeffs, len1,
                                   op2->coeffs, len2, ctx);
        _TEMPLATE(T, poly_set_length) (rop, rlen, ctx);
    }
}


#endif
