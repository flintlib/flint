/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_set_ui(qqbar_t res, ulong x)
{
    fmpz_t t;
    fmpz_init_set_ui(t, x);
    qqbar_set_fmpz(res, t);
    fmpz_clear(t);
}

