/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "partitions.h"

/* defined in flint*/
FLINT_DLL extern const unsigned int partitions_lookup[128];

/* we get log2(p(n))/2 bits with the leading term */
static int
use_exact(const fmpz_t n, slong prec)
{
    if (fmpz_size(n) >= 3)
        return 0;
    return (prec + 20.0) * (prec + 20.0) > 3.42 * fmpz_get_d(n);
}

void
arb_partitions_fmpz(arb_t res, const fmpz_t n, slong prec)
{
    if (fmpz_cmp_ui(n, 128) < 0)
    {
        arb_set_ui(res, fmpz_sgn(n) >= 0 ? partitions_lookup[*n] : 0);
        arb_set_round(res, res, prec);
    }
    else if (use_exact(n, prec))
    {
        fmpz_t t;
        fmpz_init(t);
        partitions_fmpz_fmpz(t, n, 0);
        arb_set_round_fmpz(res, t, prec);
        fmpz_clear(t);
    }
    else
    {
        mag_t err;
        mag_init(err);
        partitions_leading_fmpz(res, n, prec + 10);

        /* n >= 128; it is easy to check that the error is bounded
           by the square root of the leading approximation */
        arb_get_mag(err, res);
        mag_sqrt(err, err);
        arb_add_error_mag(res, err);

        arb_set_round(res, res, prec);
        mag_clear(err);
    }
}

void
arb_partitions_ui(arb_t res, ulong n, slong prec)
{
    if (n < 128)
    {
        arb_set_ui(res, partitions_lookup[n]);
        arb_set_round(res, res, prec);
    }
    else
    {
        fmpz_t t;
        fmpz_init_set_ui(t, n);
        arb_partitions_fmpz(res, t, prec);
        fmpz_clear(t);
    }
}

