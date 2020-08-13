/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

int 
fmpz_read(fmpz_t f)
{
    mpz_t t;
    size_t r;

    mpz_init(t);
    r = mpz_inp_str(t, stdin, 10);
    fmpz_set_mpz(f, t);
    mpz_clear(t);

    return (r > 0) ? 1 : 0;
}
