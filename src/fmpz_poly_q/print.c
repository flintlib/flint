/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/


#include "fmpz_poly_q.h"

int fmpz_poly_q_print(const fmpz_poly_q_t op)
{
    char *str;

    str = fmpz_poly_q_get_str(op);
    flint_printf("%s", str);
    flint_free(str);

    return 1;
}
