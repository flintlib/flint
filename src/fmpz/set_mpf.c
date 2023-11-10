/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "fmpz.h"

void
fmpz_set_mpf(fmpz_t f, const mpf_t x)
{
    int check;

    check = flint_mpf_fits_slong_p(x);

    if (check)
    {
        slong cx = flint_mpf_get_si(x);
        fmpz_set_si(f, cx);
    }
    else
    {
        __mpz_struct *z = _fmpz_promote(f);
        mpz_set_f(z, x);
    }
}
