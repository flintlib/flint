/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_q.h"

void fmpz_poly_q_set_si(fmpz_poly_q_t rop, slong op)
{
    fmpz_poly_set_si(rop->num, op);
    fmpz_poly_set_si(rop->den, 1);
}

