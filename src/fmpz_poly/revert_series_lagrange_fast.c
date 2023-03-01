/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011-2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "ulong_extras.h"

/* pointer to (x/Q)^i */
#define Ri(ii) (R + (n-1)*((ii)-1))

void
_fmpz_poly_revert_series_lagrange_fast(fmpz * Qinv,
    const fmpz * Q, slong Qlen, slong n)
{
    slong i, j, k, m;
    fmpz *R, *S, *T, *tmp;
    fmpz_t t;

    if (n <= 2)
    {
        _fmpz_vec_set(Qinv, Q, n);
        return;
    }

    m = n_sqrt(n);

    fmpz_init(t);
    R = _fmpz_vec_init((n - 1) * m);
    S = _fmpz_vec_init(n - 1);
    T = _fmpz_vec_init(n - 1);

    fmpz_zero(Qinv);
    fmpz_set(Qinv + 1, Q + 1);

    _fmpz_poly_inv_series(Ri(1), Q + 1, FLINT_MIN(Qlen, n) - 1, n - 1);
    for (i = 2; i <= m; i++)
        _fmpz_poly_mullow(Ri(i), Ri(i-1), n - 1, Ri(1), n - 1, n - 1);
    for (i = 2; i < m; i++)
        fmpz_divexact_ui(Qinv + i, Ri(i) + i - 1, i);

    _fmpz_vec_set(S, Ri(m), n - 1);

    for (i = m; i < n; i += m)
    {
        fmpz_divexact_ui(Qinv + i, S + i - 1, i);

        for (j = 1; j < m && i + j < n; j++)
        {
            fmpz_mul(t, S + 0, Ri(j) + i + j - 1);
            for (k = 1; k <= i + j - 1; k++)
                fmpz_addmul(t, S + k, Ri(j) + i + j - 1 - k);
            fmpz_divexact_ui(Qinv + i + j, t, i + j);
        }

        if (i + 1 < n)
        {
            _fmpz_poly_mullow(T, S, n - 1, Ri(m), n - 1, n - 1);
            tmp = S; S = T; T = tmp;
        }
    }

    fmpz_clear(t);
    _fmpz_vec_clear(R, (n - 1) * m);
    _fmpz_vec_clear(S, n - 1);
    _fmpz_vec_clear(T, n - 1);
}

void
fmpz_poly_revert_series_lagrange_fast(fmpz_poly_t Qinv,
                                        const fmpz_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    if (Qlen < 2 || !fmpz_is_zero(Q->coeffs) || !fmpz_is_pm1(Q->coeffs + 1))
    {
        flint_printf("Exception (fmpz_poly_revert_series_lagrange_fast). Input must have \n"
               "zero constant term and +1 or -1 as coefficient of x^1\n.");
        flint_abort();
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_revert_series_lagrange_fast(Qinv->coeffs, Q->coeffs, Qlen, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_revert_series_lagrange_fast(t->coeffs, Q->coeffs, Qlen, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }
    
    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);
}
