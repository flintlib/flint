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

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_factor_print(const fmpz_mod_poly_factor_t fac,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong i;

    for (i = 0; i < fac->num; i++)
    {
        fmpz_mod_poly_print(fac->poly + i, ctx);
        flint_printf(" ^ %wd\n", fac->exp[i]);
    }
}
