/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

static void
theta3_qexp(fmpz * r, slong n)
{
    slong i, j;

    _fmpz_vec_zero(r, n);

    for (i = j = 0; j < n; i++)
    {
        fmpz_set_ui(r + j, i == 0 ? 1 : 2);
        j += 1 + 2*i;
    }
}

static void
theta3_qexp_squared(fmpz * r, slong n)
{
    slong i, j, x, y;

    _fmpz_vec_zero(r, n);

    for (x = i = 0; x < n; i++)
    {
        for (y = j = 0; x + y < n; j++)
        {
            fmpz_add_ui(r + x + y, r + x + y, (x ? 2 : 1) * (y ? 2 : 1));
            y += 2 * j + 1;
        }
        x += 2 * i + 1;
    }
}

void
arith_sum_of_squares_vec(fmpz * r, ulong k, slong n)
{
    if (k == 0 || n <= 1)
    {
        _fmpz_vec_zero(r, n);
        if (n > 0)
            fmpz_set_ui(r, 1);
    }
    else if (k == 1)
    {
        theta3_qexp(r, n);
    }
    else if (k == 2)
    {
        theta3_qexp_squared(r, n);
    }
    else if (k % 2 == 0)
    {
        fmpz * t = _fmpz_vec_init(n);

        theta3_qexp_squared(t, n);
        _fmpz_poly_pow_trunc(r, t, k / 2, n);
        _fmpz_vec_clear(t, n);
    }
    else
    {
        fmpz *t, *u;
        t = _fmpz_vec_init(n);
        u = _fmpz_vec_init(n);

        theta3_qexp_squared(t, n);

        if (k == 3)
        {
            theta3_qexp(u, n);
            _fmpz_poly_mullow(r, t, n, u, n, n);
        }
        else
        {
            _fmpz_poly_pow_trunc(u, t, (k - 1) / 2, n);
            theta3_qexp(t, n);
            _fmpz_poly_mullow(r, t, n, u, n, n);
        }

        _fmpz_vec_clear(t, n);
        _fmpz_vec_clear(u, n);
    }
}
