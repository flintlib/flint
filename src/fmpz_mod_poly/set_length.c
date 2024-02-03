/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_set_length(fmpz_mod_poly_t poly, slong len)
{
    if (poly->length > len)
    {
        slong i;

        for (i = len; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
    }
    poly->length = len;
}
