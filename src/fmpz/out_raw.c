/*
    Copyright (C) 2013 Qingwen GUAN

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

size_t fmpz_out_raw( FILE *fout, const fmpz_t x ) 
{
    mpz_t v;
    size_t size;

    mpz_init( v );
    fmpz_get_mpz( v, x );
    size = mpz_out_raw( fout, v );
    mpz_clear( v );

    return size;
}


