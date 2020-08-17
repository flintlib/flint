/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"


void
_fmpz_poly_div_root(fmpz * Q, const fmpz * A, slong len, const fmpz_t c)
{
    fmpz_t r, t;
    slong i;

    if (len < 2)
        return;

    fmpz_init(r);
    fmpz_init(t);
    fmpz_set(r, A + len - 1);

    for (i = len - 2; i > 0; i--)
    {
        fmpz_mul(t, r, c);
        fmpz_add(t, t, A + i);
        fmpz_swap(Q + i, r);
        fmpz_swap(r, t);
    }

    fmpz_swap(Q, r);
    fmpz_clear(r);
    fmpz_clear(t);
}

void
fmpz_poly_div_root(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_t c)
{
    slong len = A->length;

    if (len <= 1)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (fmpz_is_zero(c))
    {
        fmpz_poly_shift_right(Q, A, 1);
        return;
    }

    fmpz_poly_fit_length(Q, len - 1);
    _fmpz_poly_div_root(Q->coeffs, A->coeffs, len, c);
    _fmpz_poly_set_length(Q, len - 1);
}
