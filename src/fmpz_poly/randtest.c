/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

void
fmpz_poly_randtest(fmpz_poly_t f, flint_rand_t state, 
                   slong len, flint_bitcnt_t bits)
{
    fmpz_poly_fit_length(f, len);
    _fmpz_vec_randtest(f->coeffs, state, len, bits);
    _fmpz_poly_set_length(f, len);
    _fmpz_poly_normalise(f);
}

void
fmpz_poly_randtest_unsigned(fmpz_poly_t f, flint_rand_t state, 
                            slong len, flint_bitcnt_t bits)
{
    fmpz_poly_fit_length(f, len);
    _fmpz_vec_randtest_unsigned(f->coeffs, state, len, bits);
    _fmpz_poly_set_length(f, len);
    _fmpz_poly_normalise(f);
}

void
fmpz_poly_randtest_not_zero(fmpz_poly_t f, flint_rand_t state, 
                            slong len, flint_bitcnt_t bits)
{
    if ((bits == 0) || (len == 0))
    {
        flint_printf("Exception (fmpz_poly_randtest_not_zero). bits or len is zero.\n");
        flint_abort();
    }

    fmpz_poly_randtest(f, state, len, bits);
    if (fmpz_poly_is_zero(f))
        fmpz_poly_set_ui(f, 1);
}
