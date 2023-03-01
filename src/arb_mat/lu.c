/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int
arb_mat_lu(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec)
{
    if (arb_mat_nrows(A) < 8 || arb_mat_ncols(A) < 8)
        return arb_mat_lu_classical(P, LU, A, prec);
    else
        return arb_mat_lu_recursive(P, LU, A, prec);
}

