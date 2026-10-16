/*
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_mulhigh(fmpz * res, const fmpz * poly1, slong len1,
                                   const fmpz * poly2, slong len2, slong start)
{
    /* the low coefficients of res are left arbitrary */
    _fmpz_poly_mulmid(res + start, poly1, len1, poly2, len2, start, len1 + len2 - 1);
}
