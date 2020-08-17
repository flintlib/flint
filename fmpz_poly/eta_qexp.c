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
_eta_one(fmpz * c, slong N)
{
    slong k, n;
    int s;

    _fmpz_vec_zero(c, N);

    /* P */
    for (k = 0, n = 0, s = 1; n < N; n += 3 * k + 1, k++, s = -s)
    {
        c[n] = s;
    }

    /* Q */
    for (k = 1, n = 2, s = -1; n < N; n += 3 * k + 2, k++, s = -s)
    {
        c[n] = s;
    }
}

void
_eta_two(fmpz * c, slong N)
{
    slong k1, n1, k2, n2;
    int s, t;

    _fmpz_vec_zero(c, N);

    /* P^2 */
    for (k1 = 0, n1 = 0; 2 * n1 < N; n1 += 3 * k1 + 1, k1++)
    {
        c[2 * n1] += 1;
    }
    for (k1 = 0, n1 = 0; n1 < N; n1 += 3 * k1 + 1, k1++)
    {
        for (k2 = k1 + 1, n2 = n1 + 3 * k1 + 1, s = -2;
                n1 + n2 < N; n2 += 3 * k2 + 1, k2++, s = -s)
        {
            c[n1 + n2] += s;
        }
    }

    /* Q^2 */
    for (k1 = 1, n1 = 2; 2 * n1 < N; n1 += 3 * k1 + 2, k1++)
    {
        c[2 * n1] += 1;
    }
    for (k1 = 1, n1 = 2; n1 < N; n1 += 3 * k1 + 2, k1++)
    {
        for (k2 = k1 + 1, n2 = n1 + 3 * k1 + 2, s = -2;
                n1 + n2 < N; n2 += 3 * k2 + 2, k2++, s = -s)
        {
            c[n1 + n2] += s;
        }
    }

    /* 2PQ */
    for (k1 = 0, n1 = 0, s = 2; n1 < N; n1 += 3 * k1 + 1, k1++, s = -s)
    {
        for (k2 = 1, n2 = 2, t = -s;
                n1 + n2 < N; n2 += 3 * k2 + 2, k2++, t = -t)
        {
            c[n1 + n2] += t;
        }
    }
}

/* R */
static void
_eta_three(fmpz * c, slong N)
{
    slong k, n;

    _fmpz_vec_zero(c, N);

    for (k = 0, n = 0; n < N; n += k + 1, k++)
        c[n] = (k % 2) ? -(2 * k + 1) : (2 * k + 1);
}

/* (P + Q) * R */
static void
_eta_four(fmpz * c, slong N)
{
    slong k1, n1, k2, n2;

    _fmpz_vec_zero(c, N);

    /* P * R */
    for (k1 = 0, n1 = 0; n1 < N; n1 += 3 * k1 + 1, k1++)
    {
        for (k2 = 0, n2 = 0; n1 + n2 < N; n2 += k2 + 1, k2++)
        {
            if ((k1 + k2) % 2)
                fmpz_sub_ui(c + n1 + n2, c + n1 + n2, 2 * k2 + 1);
            else
                fmpz_add_ui(c + n1 + n2, c + n1 + n2, 2 * k2 + 1);
        }
    }

    /* Q * R */
    for (k1 = 1, n1 = 2; n1 < N; n1 += 3 * k1 + 2, k1++)
    {
        for (k2 = 0, n2 = 0; n1 + n2 < N; n2 += k2 + 1, k2++)
        {
            if ((k1 + k2) % 2)
                fmpz_sub_ui(c + n1 + n2, c + n1 + n2, 2 * k2 + 1);
            else
                fmpz_add_ui(c + n1 + n2, c + n1 + n2, 2 * k2 + 1);
        }
    }
}

static void
_eta_six(fmpz * c, slong N)
{
    slong k1, n1, k2, n2;
    fmpz_t tmp;

    _fmpz_vec_zero(c, N);
    fmpz_init(tmp);

    /* R^2 */
    for (k1 = 0, n1 = 0; 2 * n1 < N; n1 += k1 + 1, k1++)
    {
        fmpz_set_ui(c + 2 * n1, 2 * k1 + 1);
        fmpz_mul_ui(c + 2 * n1, c + 2 * n1, 2 * k1 + 1);
    }
    for (k1 = 0, n1 = 0; n1 < N; n1 += k1 + 1, k1++)
    {
        fmpz_set_ui(tmp, 2 * (2 * k1 + 1));

        for (k2 = k1 + 1, n2 = n1 + k1 + 1; n1 + n2 < N; n2 += k2 + 1, k2++)
        {
            if ((k1 + k2) % 2)
                fmpz_submul_ui(c + n1 + n2, tmp, 2 * k2 + 1);
            else
                fmpz_addmul_ui(c + n1 + n2, tmp, 2 * k2 + 1);
        }
    }

    fmpz_clear(tmp);
}

void
_fmpz_poly_eta_qexp(fmpz * f, slong e, slong len)
{
    if (e < 0)
    {
        fmpz * t = _fmpz_vec_init(len);
        _fmpz_poly_eta_qexp(t, -e, len);
        _fmpz_poly_inv_series(f, t, len, len);
        _fmpz_vec_clear(t, len);
        return;
    }
    else if (e == 0)
    {
        _fmpz_vec_zero(f, len);
        if (len > 0)
            fmpz_set_ui(f, 1);
    }
    else if (e == 1)
    {
        _eta_one(f, len);
    }
    else if (e == 2)
    {
        _eta_two(f, len);
    }
    else if (e == 3)
    {
        _eta_three(f, len);
    }
    else if (e == 4)
    {
        _eta_four(f, len);
    }
    else if (e == 6)
    {
        _eta_six(f, len);
    }
    else
    {
        fmpz *t;

        t = _fmpz_vec_init(len);

        if (e % 6 == 0)
        {
            _eta_six(t, len);
            e /= 6;
        }
        else if (e % 4 == 0)
        {
            _eta_four(t, len);
            e /= 4;
        }
        else if (e % 3 == 0)
        {
            _eta_three(t, len);
            e /= 3;
        }
        else if (e % 2 == 0 && 0)
        {
            _eta_two(t, len);
            e /= 2;
        }
        else
        {
            _eta_one(t, len);
        }

        if (e == 2)
        {
            _fmpz_poly_sqrlow(f, t, len, len);
        }
        else if (e == 4)
        {
            _fmpz_poly_sqrlow(f, t, len, len);
            _fmpz_poly_sqrlow(t, f, len, len);
            _fmpz_vec_swap(f, t, len);
        }
        else
        {
            _fmpz_poly_pow_trunc(f, t, e, len);
        }

        _fmpz_vec_clear(t, len);
    }
}

void
fmpz_poly_eta_qexp(fmpz_poly_t f, slong e, slong n)
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
        _fmpz_poly_eta_qexp(f->coeffs, e, n);
        _fmpz_poly_set_length(f, n);
        _fmpz_poly_normalise(f);
    }
}
