/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "fq_nmod_mat.h"

int FQ_NMOD_MAT_MUL_KS_CUTOFF(slong r, slong c, const fq_nmod_ctx_t ctx)
{
    slong d2 = FLINT_MAX(0, 12 - fq_nmod_ctx_degree(ctx));

    if (2*(r + 1)*c > d2*d2)
        return 1;
    else
        return 0;
}
