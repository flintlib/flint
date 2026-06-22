/*
    Copyright (C) 2015, Elias Tsigaridas
    Copyright (C) 2016, Vincent Delecroix

    This file is part of FLINT.

    The implementation follows the function slv_poly_root_upper_bound_2exp in
    the library SLV version 0.5 by Elias Tsigaridas (in the file poly_ops.c
    lines 67-125).

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"

slong _fmpz_poly_positive_root_upper_bound_2exp_local_max(const fmpz * pol, slong len)
{
    slong b, b0, bmin;
    slong i, j, jmin = -1;
    fmpz_t tmp;

    fmpz_init(tmp);

    FLINT_ASSERT(len >= 0);

    slong * coeffs = flint_malloc((ulong)len * sizeof(slong));
    for (i = 0; i < len; i++)
        coeffs[i] = 1;

    b0 = WORD_MIN;
    int sgn = fmpz_sgn(pol + len - 1);

    for (i = len - 2; i >= 0; i--)
    {
        if (fmpz_sgn(pol + i) == 0 || fmpz_sgn(pol + i) == sgn)
            continue;

        bmin = WORD_MAX;
        for (j = i + 1; j < len; j++)
        {
            /* compare the current bound with the log (in base 2) of */
            /* (- 2^coeffs[j] * p_i / p_j) ^ (1 / (j - i))           */
            /* which equals                                          */
            /* (coeffs[j] + log(|p_i|) + log(|p_j|)) / (j - i)       */
            b = coeffs[j];

            fmpz_abs(tmp, pol + i);
            b += fmpz_clog_ui(tmp, 2);

            fmpz_abs(tmp, pol + j);
            b -= fmpz_flog_ui(tmp, 2);

            b = (b + j - i - 1) / (j - i);

            if (b < bmin)
            {
                jmin = j;
                bmin = b;
                if (bmin < b0)
                    break;
            }
        }

        b0 = FLINT_MAX(b0, bmin);

        FLINT_ASSERT(jmin >= 0);
        coeffs[jmin] ++;
    }

    fmpz_clear(tmp);
    flint_free(coeffs);

    return b0;
}

slong _fmpz_poly_positive_root_upper_bound_2exp(const fmpz * pol, slong len)
{
    return _fmpz_poly_positive_root_upper_bound_2exp_local_max(pol, len);
}

slong fmpz_poly_positive_root_upper_bound_2exp(const fmpz_poly_t pol)
{
    slong i0;

    if (fmpz_poly_is_zero(pol))
        return 0;

    i0 = 0;
    while (fmpz_is_zero(pol->coeffs + i0))
        i0++;

    return _fmpz_poly_positive_root_upper_bound_2exp_local_max(
               pol->coeffs + i0,
               fmpz_poly_length(pol) - i0);
}
