/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

int _fmpz_poly_fprint(FILE * file, const fmpz * poly, slong len)
{
    return _fmpz_vec_fprint(file, poly, len);
}

int fmpz_poly_fprint(FILE * file, const fmpz_poly_t poly)
{
    return _fmpz_vec_fprint(file, poly->coeffs, poly->length);
}

