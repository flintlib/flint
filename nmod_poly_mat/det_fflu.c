/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "perm.h"

void
nmod_poly_mat_det_fflu(nmod_poly_t det, const nmod_poly_mat_t A)
{
    slong n = nmod_poly_mat_nrows(A);

    if (n == 0)
        nmod_poly_one(det);
    else
    {
        nmod_poly_mat_t tmp;
        slong * perm;
        nmod_poly_mat_init_set(tmp, A);
        perm = _perm_init(n);

        nmod_poly_mat_fflu(tmp, det, perm, tmp, 1);
        if (_perm_parity(perm, n))
            nmod_poly_neg(det, det);

        _perm_clear(perm);
        nmod_poly_mat_clear(tmp);
    }
}
