/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "fmpz_mini.h"

void fmpz_combit(fmpz_t f, ulong i)
{
    if (!COEFF_IS_MPZ(*f))
    {
        if (i < SMALL_FMPZ_BITCOUNT_MAX)
        {
            *f ^= (WORD(1) << i);
        }
        else  /* i >= FLINT_BITS */
        {
            mpz_mock_ptr mf = _fmpz_promote_val(f);
            mpz_combit((mpz_ptr) mf, i);
            _fmpz_demote_val(f);
        }
    }
    else
    {
        mpz_mock_ptr mf = COEFF_TO_PTR(*f);
        mpz_combit((mpz_ptr) mf, i);
        _fmpz_demote_val(f);
    }
}
