/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

int
aprcl_is_prime_divisors_in_residue(const fmpz_t n, const fmpz_t s, ulong r)
{
    int result;
    ulong i;
    fmpz_t npow, nmul, fac;

    fmpz_init(fac);
    fmpz_init_set(npow, n);
    fmpz_mod(npow, npow, s); /* npow = n mod s */
    fmpz_init_set(nmul, npow);

    result = 1;
    for (i = 1; i < r; i++)
    {
        /* if find divisor then n is composite */
        if (fmpz_divisor_in_residue_class_lenstra(fac, n, npow, s) == 1)
        {
            result = 0;
            break;
        }
        /* npow = n^i mod s */
        fmpz_mul(npow, npow, nmul);
        fmpz_mod(npow, npow, s);
    }

    fmpz_clear(fac);
    fmpz_clear(npow);
    fmpz_clear(nmul);

    return result;
}

