/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

void fmpz_mod_mat_window_init(fmpz_mod_mat_t window, const fmpz_mod_mat_t mat,
                              slong r1, slong c1, slong r2, slong c2)
{
    fmpz_mat_window_init(window->mat, mat->mat, r1, c1, r2, c2);
    fmpz_init_set(window->mod, mat->mod);
}
