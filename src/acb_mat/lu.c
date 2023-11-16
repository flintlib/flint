/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
acb_mat_lu(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec)
{
    if (acb_mat_nrows(A) < 8 || acb_mat_ncols(A) < 8)
        return acb_mat_lu_classical(P, LU, A, prec);
    else
        return acb_mat_lu_recursive(P, LU, A, prec);
}

