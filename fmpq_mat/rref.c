/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

slong
fmpq_mat_rref(fmpq_mat_t B, const fmpq_mat_t A)
{
    if (A->r <= 2 || A->c <= 2)
        return fmpq_mat_rref_classical(B, A);
    else
        return fmpq_mat_rref_fraction_free(B, A);
}
