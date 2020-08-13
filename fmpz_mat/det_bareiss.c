/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "perm.h"

void
_fmpz_mat_det_bareiss(fmpz_t det, fmpz_mat_t tmp)
{
    slong *perm, n = fmpz_mat_nrows(tmp);
    perm = _perm_init(n);

    fmpz_mat_fflu(tmp, det, perm, tmp, 1);

    if (_perm_parity(perm, n) == 1)
        fmpz_neg(det, det);

    _perm_clear(perm);
}

void
fmpz_mat_det_bareiss(fmpz_t det, const fmpz_mat_t A)
{
    fmpz_mat_t tmp;

    if (A->r < 1)
    {
        fmpz_one(det);
        return;
    }

    fmpz_mat_init_set(tmp, A);
    _fmpz_mat_det_bareiss(det, tmp);
    fmpz_mat_clear(tmp);
}
