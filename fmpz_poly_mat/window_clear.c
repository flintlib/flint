/*
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2015 Elena Sergeichave

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_window_clear(fmpz_poly_mat_t window)
{
    if (window->r != 0)
        flint_free(window->rows);
}
