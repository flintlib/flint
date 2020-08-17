/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "padic_mat.h"

int padic_mat_equal(const padic_mat_t A, const padic_mat_t B)
{
    if (A->val == B->val)
    {
        return fmpz_mat_equal(padic_mat(A), padic_mat(B));
    }
    else
    {
        return 0;
    }
}

