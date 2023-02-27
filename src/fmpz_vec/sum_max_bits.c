/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"

/*
    set sumabs = bit_count(sum of absolute values of coefficients);
    set maxabs = bit_count(max of absolute values of coefficients);
*/
void _fmpz_vec_sum_max_bits(slong * sumabs, slong * maxabs,
                                             const fmpz * coeffs, slong length)
{
    slong j;
    ulong hi = 0, lo = 0;

    maxabs[0] = 0;

    for (j = 0; j < length && fmpz_fits_si(coeffs + j); j++)
    {
        slong c = fmpz_get_si(coeffs + j);
        ulong uc = (ulong) FLINT_ABS(c);
        add_ssaaaa(hi, lo, hi, lo, UWORD(0), uc);
        maxabs[0] = FLINT_MAX(maxabs[0], FLINT_BIT_COUNT(uc));
    }

    if (j == length)
    {
        /* no large coeffs encountered */
        if (hi != 0)
            sumabs[0] = FLINT_BIT_COUNT(hi) + FLINT_BITS;
        else
            sumabs[0] = FLINT_BIT_COUNT(lo);
    } else
    {
        /* restart with multiprecision routine */
        fmpz_t sum;
        fmpz_init(sum);

        for (j = 0; j < length; j++)
        {
            slong this_size = fmpz_sizeinbase(coeffs + j, 2);

            maxabs[0] = FLINT_MAX(maxabs[0], this_size);

            if (fmpz_sgn(coeffs + j) < 0)
                fmpz_sub(sum, sum, coeffs + j);
            else
                fmpz_add(sum, sum, coeffs + j);
        }

        sumabs[0] = fmpz_sizeinbase(sum, 2);
        fmpz_clear(sum);
   }
}
