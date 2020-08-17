/*
    Copyright (C) 2017 Tommy Hofmann 

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define FMPZ_FACTOR_INLINES_C

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void fmpz_factor_get_fmpz(fmpz_t z, const fmpz_factor_t factor, slong i)
{
    fmpz_set(z, factor->p + i);
}

