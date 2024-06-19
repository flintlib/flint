/*
    Copyright (C) 2012 Andres Goens
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

void
_TEMPLATE(T, poly_rem)(TEMPLATE(T, struct) *R,
                       const TEMPLATE(T, struct) *A, slong lenA,
                       const TEMPLATE(T, struct) *B, slong lenB,
                       const TEMPLATE(T, t) invB,
                       const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) *Q = _TEMPLATE(T, vec_init)(lenA - lenB + 1, ctx);

    if (lenA < lenB)
    {
        _TEMPLATE(T, vec_set)(R, A, lenA, ctx);
        _TEMPLATE(T, vec_zero)(R + lenA, lenB - 1 - lenA, ctx);
    }
    else
    {
        TEMPLATE(T, struct) *T = _TEMPLATE(T, vec_init)(lenA, ctx);
       _TEMPLATE(T, poly_divrem)(Q, T, A, lenA, B, lenB, invB, ctx);
       _TEMPLATE(T, vec_set)(R, T, lenB - 1, ctx);
       _TEMPLATE(T, vec_clear)(T, lenA, ctx);
    }

    _TEMPLATE(T, vec_clear)(Q, lenA - lenB + 1, ctx);
}

void
TEMPLATE(T, poly_rem)(TEMPLATE(T, poly_t) R,
                      const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B,
                      const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) Q;
    TEMPLATE(T, poly_init)(Q, ctx);
    TEMPLATE(T, poly_divrem)(Q, R, A, B, ctx);
    TEMPLATE(T, poly_clear)(Q, ctx);
}

#endif
