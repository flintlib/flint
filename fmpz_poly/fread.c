/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz.h"
#include "fmpz_poly_mini.h"
#ifdef LONGSLONG
# define flint_mpz_get_si mpz_get_si
#else
# include "gmpcompat.h"
#endif

int fmpz_poly_fread(FILE * file, fmpz_poly_t poly)
{
    int r;
    slong i, len;
    mpz_t t;

    mpz_init(t);
    r = mpz_inp_str(t, file, 10);
    if (r == 0)
    {
        mpz_clear(t);
        return 0;
    }
    if (!mpz_fits_slong_p(t))
        flint_throw(FLINT_ERROR, "Length does not fit into a slong in fmpz_poly_fread\n");

    len = flint_mpz_get_si(t);
    mpz_clear(t);

    fmpz_poly_fit_length(poly, len);

    for (i = 0; i < len; i++)
    {
        r = fmpz_fread(file, poly->coeffs + i);
        if (r <= 0)
            return r;
    }

    _fmpz_poly_set_length(poly, len);
    _fmpz_poly_normalise(poly);

    return 1;
}
