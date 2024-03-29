/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void TEMPLATE(T, gcdinv)(TEMPLATE(T, t) rop, TEMPLATE(T, t) inv,
                    const TEMPLATE(T, t) op,
                    const TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, is_invertible)(const TEMPLATE(T, t) op,
                           const TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, is_invertible_f)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                             const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, div)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op1,
                 const TEMPLATE(T, t) op2, const TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, multiplicative_order)(fmpz * ord, const TEMPLATE(T, t) op,
                             const TEMPLATE(T, ctx_t) ctx);

#ifdef B
void TEMPLATE4(T, get, B, mat)(TEMPLATE(B, mat_t) col,
                                         const TEMPLATE(T, t) a,
                                         const TEMPLATE(T, ctx_t) ctx);
void TEMPLATE4(T, set, B, mat)(TEMPLATE(T, t) a,
                                         const TEMPLATE(B, mat_t) col,
                                         const TEMPLATE(T, ctx_t) ctx);
#endif


FQ_TEMPLATES_INLINE
int TEMPLATE(T, is_primitive)(const TEMPLATE(T, t) op, const TEMPLATE(T, ctx_t) ctx)
{
    int ret;
    ret = TEMPLATE(T, multiplicative_order)(NULL, op, ctx) == 1;
    return ret;
}

#endif
