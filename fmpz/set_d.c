/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz-impl.h"

void
fmpz_set_d(fmpz_t f, double c)
{
    if (c >= DOUBLE_MIN && c <= DOUBLE_MAX)
    {
        _fmpz_demote(f);
        /* guaranteed to fit, since c gets truncated */
        *f = (slong) c;
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(f);
        mpz_set_d(z, c);
        _fmpz_demote_val(f);
    }
}
