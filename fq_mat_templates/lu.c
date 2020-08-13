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

slong
TEMPLATE(T, mat_lu) (slong * P, TEMPLATE(T, mat_t) A, int rank_check,
                     const TEMPLATE(T, ctx_t) ctx)
{
    return TEMPLATE(T, mat_lu_recursive) (P, A, rank_check, ctx);
}


#endif
