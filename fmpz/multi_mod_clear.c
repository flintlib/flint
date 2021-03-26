/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"


void fmpz_multi_mod_clear(fmpz_multi_mod_t P)
{
    slong i;

    for (i = 0; i < P->alloc; i++)
    {
        fmpz_clear(P->prog[i].modulus);
        fmpz_clear(P->moduli + i);
    }

    flint_free(P->prog);
    flint_free(P->moduli);
}

