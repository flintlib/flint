/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_fit_length(fmpz_mod_poly_t f, slong len)
{
    if (f->alloc < len)
    {
        slong alloc = FLINT_MAX(2*f->alloc, len);
        f->coeffs = FLINT_ARRAY_REALLOC(f->coeffs, alloc, fmpz);
        flint_mpn_zero((mp_ptr) (f->coeffs + f->alloc), alloc - f->alloc);
        f->alloc = alloc;
    }
}

