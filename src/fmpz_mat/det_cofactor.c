/*
    Copyright (C) 2010,2011,2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

#define E(i,j) fmpz_mat_entry(A, i, j)

static void
_fmpz_mat_det_cofactor_2x2(fmpz_t det, const fmpz_mat_t A)
{
    fmpz_fmms(det, E(0,0), E(1,1), E(0,1), E(1,0));
}

static void
_fmpz_mat_det_cofactor_3x3(fmpz_t det, const fmpz_mat_t A)
{
    fmpz_t a;
    fmpz_init(a);

    fmpz_fmms(a, E(1,0), E(2,1), E(1,1), E(2,0));
    fmpz_mul(det, a, E(0,2));
    fmpz_fmms(a, E(1,2), E(2,0), E(1,0), E(2,2));
    fmpz_addmul(det, a, E(0,1));
    fmpz_fmms(a, E(1,1), E(2,2), E(1,2), E(2,1));
    fmpz_addmul(det, a, E(0,0));

    fmpz_clear(a);
}

static void
_fmpz_mat_det_cofactor_4x4(fmpz_t det, const fmpz_mat_t A)
{
    fmpz_t a, b;
    fmpz_init(a);
    fmpz_init(b);

    fmpz_fmms(a, E(0,3), E(1,2), E(0,2), E(1,3));
    fmpz_fmms(b, E(2,1), E(3,0), E(2,0), E(3,1));
    fmpz_mul(det, a, b);

    fmpz_fmms(a, E(0,1), E(1,3), E(0,3), E(1,1));
    fmpz_fmms(b, E(2,2), E(3,0), E(2,0), E(3,2));
    fmpz_addmul(det, a, b);

    fmpz_fmms(a, E(0,2), E(1,1), E(0,1), E(1,2));
    fmpz_fmms(b, E(2,3), E(3,0), E(2,0), E(3,3));
    fmpz_addmul(det, a, b);

    fmpz_fmms(a, E(0,3), E(1,0), E(0,0), E(1,3));
    fmpz_fmms(b, E(2,2), E(3,1), E(2,1), E(3,2));
    fmpz_addmul(det, a, b);

    fmpz_fmms(a, E(0,0), E(1,2), E(0,2), E(1,0));
    fmpz_fmms(b, E(2,3), E(3,1), E(2,1), E(3,3));
    fmpz_addmul(det, a, b);

    fmpz_fmms(a, E(0,1), E(1,0), E(0,0), E(1,1));
    fmpz_fmms(b, E(2,3), E(3,2), E(2,2), E(3,3));
    fmpz_addmul(det, a, b);

    fmpz_clear(a);
    fmpz_clear(b);
}

void
fmpz_mat_det_cofactor(fmpz_t det, const fmpz_mat_t A)
{
    switch (fmpz_mat_nrows(A))
    {
        case 0:  fmpz_one(det);                          break;
        case 1:  fmpz_set(det, fmpz_mat_entry(A, 0, 0)); break;
        case 2:  _fmpz_mat_det_cofactor_2x2(det, A);     break;
        case 3:  _fmpz_mat_det_cofactor_3x3(det, A);     break;
        case 4:  _fmpz_mat_det_cofactor_4x4(det, A);     break;
        default: flint_throw(FLINT_ERROR, "Exception (fmpz_mat_det_cofactor). dim > 4 not implemented.");
    }
}
