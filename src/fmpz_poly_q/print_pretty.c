/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>

#include "fmpz_poly_q.h"

int fmpz_poly_q_print_pretty(const fmpz_poly_q_t op, const char *x)
{
    char *str;

    str = fmpz_poly_q_get_str_pretty(op, x);
    flint_printf("%s", str);
    flint_free(str);

    return 1;
}
