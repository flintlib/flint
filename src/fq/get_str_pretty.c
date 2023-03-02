/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

char *
fq_get_str_pretty(const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_get_str_pretty(op, ctx->var);
}
