/*
    Copyright (C) 2019 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

slong fmpz_mod_mat_howell_form(fmpz_mod_mat_t mat)
{
    return fmpz_mat_howell_form_mod(mat->mat, mat->mod);
}
