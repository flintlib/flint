/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"

void _fmpz_mod_vec_dot_rev(fmpz_t r, const fmpz * a,
		           const fmpz * b, slong len, const fmpz_mod_ctx_t ctx)
{
    slong i;
    
    fmpz_zero(r);

    for (i = 0; i < len; i++)
        fmpz_addmul(r, a + i, b + len - i - 1);

    fmpz_mod_set_fmpz(r, r, ctx);
}
