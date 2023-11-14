/*
    Copyright 2021, 2022 William Hart
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

static inline
int _TEMPLATE(T, poly_sqrt_2)(TEMPLATE(T, struct) * s,
                   const TEMPLATE(T, struct) * p, slong len, TEMPLATE(T, ctx_t) ctx)
{
   slong i;

   for (i = 1; i < len; i += 2)
       if (!TEMPLATE(T, is_zero)(p + i, ctx))
           return 0;

   for (i = 0; i < len; i += 2)
       TEMPLATE(T, sqrt)(s + i / 2,  p + i, ctx);

   return 1;
}

int
_TEMPLATE(T, poly_sqrt)(TEMPLATE(T, struct) * s, const TEMPLATE(T, struct) * p, slong len, TEMPLATE(T, ctx_t) ctx)
{
    slong slen, i;
    int result;
    TEMPLATE(T, struct) * t;
    TEMPLATE(T, t) c, d;

    if (len % 2 == 0)
        return len == 0;

    if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) == 0)
        return _TEMPLATE(T, poly_sqrt_2)(s, p, len, ctx);

    /* valuation must be even, and then can be reduced to 0 */
    while (TEMPLATE(T, is_zero)(p + 0, ctx))
    {
        if (!TEMPLATE(T, is_zero)(p + 1, ctx))
            return 0;

        TEMPLATE(T, zero)(s + 0, ctx);
        p += 2;
        len -= 2;
        s++;
    }

    TEMPLATE(T, init)(c, ctx);
    TEMPLATE(T, init)(d, ctx);

    TEMPLATE(T, set)(d, p + 0, ctx);
    TEMPLATE(T, set)(c, d, ctx);

    if (!TEMPLATE(T, is_one)(c, ctx))
    {
        if (!TEMPLATE(T, sqrt)(c, c, ctx))
        {
           result = 0;
           goto cleanup;
        }
    }

    if (len == 1)
    {
        TEMPLATE(T, set)(s + 0, c, ctx);

        result = 1;
        goto cleanup;
    }

    slen = len / 2 + 1;

    t = _TEMPLATE(T, vec_init)(len, ctx);

    if (TEMPLATE(T, is_one)(c, ctx))
        _TEMPLATE(T, poly_sqrt_series)(s, p, slen, ctx);
    else
    {
        TEMPLATE(T, inv)(d, d, ctx);

        _TEMPLATE3(T, vec_scalar_mul, T)(t, p, slen, d, ctx);
        _TEMPLATE(T, poly_sqrt_series)(s, t, slen, ctx);
    }

    if (!TEMPLATE(T, is_one)(c, ctx))
        _TEMPLATE3(T, vec_scalar_mul, T)(s, s, slen, c, ctx);

    _TEMPLATE(T, poly_mulhigh)(t, s, slen, s, slen, slen, ctx);

    for (i = 0; i < slen; i++)
        TEMPLATE(T, zero)(t + i, ctx);

    result = _TEMPLATE(T, vec_equal)(t + slen, p + slen, len - slen, ctx);

    _TEMPLATE(T, vec_clear)(t, len, ctx);

cleanup:

    TEMPLATE(T, clear)(c, ctx);
    TEMPLATE(T, clear)(d, ctx);

    return result;
}

int
TEMPLATE(T, poly_sqrt)(TEMPLATE(T, poly_t) b, const TEMPLATE(T, poly_t) a, TEMPLATE(T, ctx_t) ctx)
{
    slong blen, len = a->length;
    int result;

    if (len % 2 == 0)
    {
        TEMPLATE(T, poly_zero)(b, ctx);

        return len == 0;
    }

    if (b == a)
    {
        TEMPLATE(T, poly_t) tmp;

        TEMPLATE(T, poly_init)(tmp, ctx);

        result = TEMPLATE(T, poly_sqrt)(tmp, a, ctx);
        TEMPLATE(T, poly_swap)(b, tmp, ctx);

        TEMPLATE(T, poly_clear)(tmp, ctx);

        return result;
    }

    blen = len / 2 + 1;
    TEMPLATE(T, poly_fit_length)(b, blen, ctx);

    result = _TEMPLATE(T, poly_sqrt)(b->coeffs, a->coeffs, len, ctx);

    if (!result)
        blen = 0;

    _TEMPLATE(T, poly_set_length)(b, blen, ctx);
    _TEMPLATE(T, poly_normalise)(b, ctx);

    return result;
}

#endif

