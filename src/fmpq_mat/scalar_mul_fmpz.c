/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"

void fmpq_mat_scalar_mul_fmpz(fmpq_mat_t rop,
                              const fmpq_mat_t op, const fmpz_t x)
{
    slong i, j;

    for (i = 0; i < op->r; i++)
        for (j = 0; j < op->c; j++)
            fmpq_mul_fmpz(fmpq_mat_entry(rop, i, j),
                          fmpq_mat_entry(op, i, j), x);
}
