/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"

/* todo: efficient implementation */
char *
_flint_mpn_get_str(mp_srcptr x, mp_size_t n)
{
    fmpz_t t;
    char * s;
    fmpz_init(t);
    fmpz_set_ui_array(t, x, n);
    s = fmpz_get_str(NULL, 10, t);
    fmpz_clear(t);
    return s;
}