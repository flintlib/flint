/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "long_extras.h"

/*
    If |f|_infty <= A and degs hold the degrees of f, set B to a bound on
    |g|_infty for any divisor g of f. Return 1 for success, 0 for failure.
*/
int fmpz_mpoly_factor_bound_si(fmpz_t B, const fmpz_t A,
                                               const slong * degs, slong nvars)
{
    slong i, n = 0;
    fmpz_t t;

    FLINT_ASSERT(nvars > 0);

    fmpz_init_set_ui(t, 1);
    for (i = 1; i < nvars; i++)
    {
        if (degs[i] < 0)
        {
            fmpz_clear(t);
            fmpz_zero(B);
            return 1;
        }
        fmpz_mul_ui(t, t, degs[i] + 1);
        if (z_add_checked(&n, n, degs[i]))
        {
            fmpz_clear(t);
            return 0;
        }
    }
    fmpz_cdiv_q_2exp(t, t, nvars);
    fmpz_sqrt(t, t);
    fmpz_add_ui(t, t, 1);
    fmpz_mul(B, A, t);
    fmpz_mul_2exp(B, B, n);
    fmpz_abs(B, B);
    fmpz_clear(t);
    return 1;
}
