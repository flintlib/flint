/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "fmpz_mini.h"
#ifdef LONGSLONG
# define flint_mpz_bin_uiui mpz_bin_uiui
#else
# include "gmpcompat.h"
#endif

/* TODO: speedup for small n,k */
void fmpz_bin_uiui(fmpz_t res, ulong n, ulong k)
{
    mpz_mock_ptr t = _fmpz_promote(res);
    flint_mpz_bin_uiui((mpz_ptr) t, n, k);
    _fmpz_demote_val(res);
}
