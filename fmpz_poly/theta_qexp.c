/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

static void
theta_one(fmpz * r, slong n)
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
theta_two(fmpz * r, slong n)
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
_fmpz_poly_theta_qexp(fmpz * f, slong k, slong n)
{
    if (k < 0)
    {
        fmpz * t = _fmpz_vec_init(n);
        _fmpz_poly_theta_qexp(t, -k, n);
        _fmpz_poly_inv_series(f, t, n, n);
        _fmpz_vec_clear(t, n);
        return;
    }
    else if (k == 0)
    {
        _fmpz_vec_zero(f, n);
        if (n > 0)
            fmpz_set_ui(f, 1);
    }
    else if (k == 1)
    {
        theta_one(f, n);
    }
    else if (k == 2)
    {
        theta_two(f, n);
    }
    else if (k % 2 == 0)
    {
        fmpz * t = _fmpz_vec_init(n);
        theta_two(t, n);
        _fmpz_poly_pow_trunc(f, t, k / 2, n);
        _fmpz_vec_clear(t, n);
    }
    else
    {
        fmpz *t, *u;
        t = _fmpz_vec_init(n);
        u = _fmpz_vec_init(n);

        theta_two(t, n);

        if (k == 3)
        {
            theta_one(u, n);
            _fmpz_poly_mullow(f, t, n, u, n, n);
        }
        else
        {
            _fmpz_poly_pow_trunc(u, t, (k - 1) / 2, n);
            theta_one(t, n);
            _fmpz_poly_mullow(f, t, n, u, n, n);
        }

        _fmpz_vec_clear(t, n);
        _fmpz_vec_clear(u, n);
    }
}

void
fmpz_poly_theta_qexp(fmpz_poly_t f, slong e, slong n)
{
    if (n < 1)
    {
        fmpz_poly_zero(f);
    }
    else if (e == 0 || n == 1)
    {
        fmpz_poly_one(f);
    }
    else
    {
        fmpz_poly_fit_length(f, n);
        _fmpz_poly_theta_qexp(f->coeffs, e, n);
        _fmpz_poly_set_length(f, n);
        _fmpz_poly_normalise(f);
    }
}
