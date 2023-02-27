/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_q.h"

void fmpz_poly_q_clear(fmpz_poly_q_t rop)
{
    if (rop->num != NULL)
    {
        fmpz_poly_clear(rop->num);
        flint_free(rop->num);
        rop->num = NULL;
    }
    if (rop->den != NULL)
    {
        fmpz_poly_clear(rop->den);
        flint_free(rop->den);
        rop->den = NULL;
    }
}

