/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "ulong_extras.h"
#include "fmpz_mini.h"

void fmpz_sqrt(fmpz_t f, const fmpz_t g)
{
    if (fmpz_sgn(g) < 0)
        flint_throw(FLINT_ERROR, "g is negative in fmpz_sqrt\n");
    
    if (!COEFF_IS_MPZ(*g))
        fmpz_set_ui(f, n_sqrt(*g));
    else
    {
        mpz_mock_ptr mf = _fmpz_promote(f);
        mpz_sqrt((mpz_ptr) mf, (mpz_ptr) COEFF_TO_PTR(*g));
        _fmpz_demote_val(f);
    }
}
