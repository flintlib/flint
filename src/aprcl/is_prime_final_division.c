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
aprcl_is_prime_final_division(const fmpz_t n, const fmpz_t s, ulong r)
{
    int result = 1;
    ulong i;
    fmpz_t npow, nmul, rem;

    fmpz_init(rem);
    fmpz_init_set(npow, n);
    fmpz_mod(npow, npow, s); /* npow = n mod s */
    fmpz_init_set(nmul, npow);

    for (i = 1; i <= r; i++)
    {
        if (fmpz_is_one(npow))
            break;

        fmpz_mod(rem, n, npow);

        /* if npow | n */
        if (fmpz_is_zero(rem))
        {
            /* if npow != n and npow != 1 */
            if (!fmpz_equal(n, npow) && !fmpz_is_one(npow))
            {
                /* npow | n, so n is composite */
                result = 0;
                break;
            }
        }

        /* npow = n^i mod s */
        fmpz_mul(npow, npow, nmul);
        fmpz_mod(npow, npow, s);
    }

    fmpz_clear(npow);
    fmpz_clear(nmul);
    fmpz_clear(rem);

    return result;
}
