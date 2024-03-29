/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_q.h"

void fmpz_poly_q_init(fmpz_poly_q_t rop)
{
    rop->num = flint_malloc(sizeof(fmpz_poly_struct));
    rop->den = flint_malloc(sizeof(fmpz_poly_struct));
    fmpz_poly_init(rop->num);
    fmpz_poly_init(rop->den);
    fmpz_poly_set_si(rop->den, 1);
}
