/*
    Copyright (C) 2010 Fredrik Johansson
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
TEMPLATE(T, mat_swap) (TEMPLATE(T, mat_t) mat1, TEMPLATE(T, mat_t) mat2,
                       const TEMPLATE(T, ctx_t) ctx)
{
    FLINT_SWAP(TEMPLATE(T, mat_struct), *mat1, *mat2);
}


#endif
