/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

/*
    The solutions for x in a*x + b*y = c && 0 <= x < b  will be precisely
        {xstart + xstride*i | 0 <= i < xlength }
*/
void fmpz_divides_mod_list(fmpz_t xstart, fmpz_t xstride, fmpz_t xlength,
                                const fmpz_t c, const fmpz_t a, const fmpz_t b)
{
    fmpz_t d, x, q, r, bbar;

    fmpz_init(d);
    fmpz_init(x);
    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(bbar);

    FLINT_ASSERT(fmpz_sgn(b) > 0);

    fmpz_mod(r, a, b);
    fmpz_gcdinv(d, x, r, b);
    fmpz_fdiv_qr(q, r, c, d);
    if (fmpz_is_zero(r))
    {
        fmpz_divexact(bbar, b, d);
        fmpz_mul(x, x, q);
        fmpz_fdiv_q(r, x, bbar);
        fmpz_submul(x, bbar, r);
        fmpz_swap(xlength, d);
        fmpz_swap(xstride, bbar);
        fmpz_swap(xstart, x);
    }
    else
    {
        fmpz_zero(xstart);
        fmpz_zero(xstride);
        fmpz_zero(xlength);
    }

    fmpz_clear(d);
    fmpz_clear(x);
    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(bbar);
}
