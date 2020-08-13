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

int padic_mat_is_zero(const padic_mat_t A)
{
    if (A->val == 0)
    {
        return fmpz_mat_is_zero(padic_mat(A));
    }
    else
    {
        return 0;
    }
}

