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
fmpz_mat_scalar_submul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
{
    if (c > 0)
        fmpz_mat_scalar_submul_ui(B, A, c);
    else
        fmpz_mat_scalar_addmul_ui(B, A, -c);
}
