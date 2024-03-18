/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_vec.h"
#include "fmpz_mod_vec.h"

void _fmpz_mod_vec_dot(fmpz_t d, const fmpz * a,
                   const fmpz * b, slong len, const fmpz_mod_ctx_t ctx)
{
    _fmpz_vec_dot(d, a, b, len);
    fmpz_mod_set_fmpz(d, d, ctx);
}

void _fmpz_mod_vec_dot_rev(fmpz_t r, const fmpz * a,
		           const fmpz * b, slong len, const fmpz_mod_ctx_t ctx)
{
    _fmpz_vec_dot_general(r, NULL, 0, a, b, 1, len);
    fmpz_mod_set_fmpz(r, r, ctx);
}
