/*
    Copyright (C) 2008, 2009, 2016 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

slong
fmpz_mpoly_get_coeff_si(const fmpz_mpoly_t poly,
                                           slong n, const fmpz_mpoly_ctx_t ctx)
{
    return (n < poly->length) ? fmpz_get_si(poly->coeffs + n) : WORD(0);
}
