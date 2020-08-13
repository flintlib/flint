/*
    Copyright (C) 2011 Fredrik Johansson
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
TEMPLATE(T, mat_submul) (TEMPLATE(T, mat_t) D,
                         const TEMPLATE(T, mat_t) C,
                         const TEMPLATE(T, mat_t) A,
                         const TEMPLATE(T, mat_t) B,
                         const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, mat_t) tmp;
    TEMPLATE(T, mat_init) (tmp, A->r, B->c, ctx);
    TEMPLATE(T, mat_mul) (tmp, A, B, ctx);
    TEMPLATE(T, mat_sub) (D, C, tmp, ctx);
    TEMPLATE(T, mat_clear) (tmp, ctx);
}


#endif
