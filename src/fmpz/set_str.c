/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int fmpz_set_str(fmpz_t f, const char * str, int b)
{
    int ans;
    mpz_t copy;

    ans = mpz_init_set_str(copy, (char *) str, b);
    if (ans == 0)
        fmpz_set_mpz(f, copy);
    mpz_clear(copy);
    return ans;
}

