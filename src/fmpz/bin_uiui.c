/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"

/* TODO: speedup for small n,k */
void fmpz_bin_uiui(fmpz_t res, ulong n, ulong k)
{
    mpz_ptr t = _fmpz_promote(res);
    flint_mpz_bin_uiui(t, n, k);
    _fmpz_demote_val(res);
}
