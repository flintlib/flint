/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
fmpz_is_probabprime(const fmpz_t p)
{
    fmpz c = *p;
    
    if (fmpz_sgn(p) <= 0)
       return 0;

    if (!COEFF_IS_MPZ(c))
       return n_is_probabprime(c);
    else
    {
       int ret;
       
       ret = (mpz_probab_prime_p(COEFF_TO_PTR(c), 25) != 0);
       return ret;
    }
}
