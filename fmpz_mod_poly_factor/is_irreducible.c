/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
fmpz_mod_poly_is_irreducible(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    if (fmpz_mod_poly_length(f, ctx) > 2)
    {
        return fmpz_mod_poly_is_irreducible_ddf(f, ctx);
    }

    return 1;
}
