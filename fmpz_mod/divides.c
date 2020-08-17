/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

int fmpz_mod_divides(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    int success;
    fmpz_t g, x, y, q;

    if (fmpz_is_zero(c))
    {
        if (fmpz_is_zero(b))
        {
            fmpz_zero(a);
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else if (fmpz_is_zero(b))
    {
        fmpz_zero(a);
        return 1;
    }

    /* b and c both not zero now */

    fmpz_init(g);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(q);

    /* solve g = c*x + n*y where g = gcd(c, n) */
    fmpz_xgcd(g, x, y, c, ctx->n);

    fmpz_fdiv_qr(q, y, b, g);
    success = fmpz_is_zero(y);
    if (success)
    {
        fmpz_mul(a, q, x);
        fmpz_mod(a, a, ctx->n);        
    }

    fmpz_clear(g);
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(q);

    return success;
}
