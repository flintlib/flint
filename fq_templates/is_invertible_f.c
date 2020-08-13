/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int
TEMPLATE(T, is_invertible_f)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                             const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) inv;
    TEMPLATE(T, init)(inv, ctx);
    TEMPLATE(T, gcdinv)(rop, inv, op, ctx);
    TEMPLATE(T, clear)(inv, ctx);
    return TEMPLATE(T, is_one)(rop, ctx);
}

#endif
