/*
    Copyright (C) 2011, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "padic_mat.h"

void padic_mat_get_entry_padic(padic_t rop, 
                               const padic_mat_t op, slong i, slong j, 
                               const padic_ctx_t ctx)
{
    fmpz_set(padic_unit(rop), padic_mat_entry(op, i, j));
    padic_val(rop) = padic_mat_val(op);
    padic_reduce(rop, ctx);
}

