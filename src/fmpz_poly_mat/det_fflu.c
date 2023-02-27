/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"
#include "perm.h"

void
fmpz_poly_mat_det_fflu(fmpz_poly_t det, const fmpz_poly_mat_t A)
{
    slong n = fmpz_poly_mat_nrows(A);

    if (n == 0)
        fmpz_poly_one(det);
    else
    {
        fmpz_poly_mat_t tmp;
        slong * perm;
        fmpz_poly_mat_init_set(tmp, A);
        perm = _perm_init(n);

        fmpz_poly_mat_fflu(tmp, det, perm, tmp, 1);
        if (_perm_parity(perm, n))
            fmpz_poly_neg(det, det);

        _perm_clear(perm);
        fmpz_poly_mat_clear(tmp);
    }
}
