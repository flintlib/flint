/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default.h"
#include "fq_default_mat.h"

void fq_default_mat_entry_set_fmpz(fq_default_mat_t mat, slong i, slong j, const fmpz_t x, const fq_default_ctx_t ctx)
{
   fq_default_t c;
   fq_default_init(c, ctx);
   fq_default_set_fmpz(c, x, ctx);
   fq_default_mat_entry_set(mat, i, j, c, ctx);
   fq_default_clear(c, ctx);
}
