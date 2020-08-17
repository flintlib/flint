/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_scalar_addmul_nmod_mat_ui(fmpz_mat_t B,
                        const nmod_mat_t A, ulong c)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_ui(t, c);
    fmpz_mat_scalar_addmul_nmod_mat_fmpz(B, A, t);
    fmpz_clear(t);
}
