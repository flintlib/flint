/*
    Copyright (C) 2013 Qingwen GUAN

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

size_t fmpz_inp_raw( fmpz_t x, FILE *fin ) 
{
    mpz_t v;
    size_t size;

    mpz_init( v );
    size = mpz_inp_raw( v, fin );
    fmpz_set_mpz( x, v );
    mpz_clear( v );

    return size;
}

