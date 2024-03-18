/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"
#include "fq_mat.h"

int FQ_MAT_MUL_KS_CUTOFF(slong r, slong c, const fq_ctx_t ctx)
{
    if (5 * FLINT_MIN(r, c) > 8 * fq_ctx_degree(ctx) + 29)
        return 1;
    else
        return 0;
}
