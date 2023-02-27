/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
fmpz_factor_si(fmpz_factor_t factor, slong n)
{
    _fmpz_factor_set_length(factor, 0);

    if (n < 0)
    {
        _fmpz_factor_extend_factor_ui(factor, -n);
        factor->sign = -1;
        return;
    }
    else
    {
        factor->sign = 1;
        _fmpz_factor_extend_factor_ui(factor, n);
    }
}
