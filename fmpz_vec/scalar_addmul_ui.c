/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_vec_scalar_addmul_ui(fmpz * vec1, const fmpz * vec2, slong len2, ulong c)
{
    ulong ix;

    for (ix = 0; ix < len2; ix++)
        fmpz_addmul_ui(vec1 + ix, vec2 + ix, c);
}
