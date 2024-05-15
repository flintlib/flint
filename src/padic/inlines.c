/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define PADIC_INLINES_C

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "padic.h"

void padic_get_unit(fmpz_t f, padic_t p)
{
   fmpz_set(f, padic_unit(p));
}
