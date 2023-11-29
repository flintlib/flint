/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2019 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_div_series(fmpz * Q, const fmpz * A, slong Alen,
    const fmpz * B, slong Blen, slong n)
{
    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);

    if (n < 72 || Blen < 72 || Alen == 1)
       _fmpz_poly_div_series_basecase(Q, A, Alen, B, Blen, n);
    else if (fmpz_is_pm1(B + 0))
    {
        fmpz * Binv = _fmpz_vec_init(n);

        _fmpz_poly_inv_series(Binv, B, Blen, n);
        _fmpz_poly_mullow(Q, Binv, n, A, Alen, n);

        _fmpz_vec_clear(Binv, n);
    } else
        _fmpz_poly_div_series_divconquer(Q, A, Alen, B, Blen, n);
}

void fmpz_poly_div_series(fmpz_poly_t Q, const fmpz_poly_t A,
                                         const fmpz_poly_t B, slong n)
{
    slong Alen = FLINT_MIN(A->length, n);
    slong Blen = FLINT_MIN(B->length, n);

    if (Blen == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_div_series). Division by zero.\n");
    }

    if (Alen == 0)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_div_series(t->coeffs, A->coeffs, Alen, B->coeffs, Blen, n);
        fmpz_poly_swap(Q, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(Q, n);
        _fmpz_poly_div_series(Q->coeffs, A->coeffs, Alen, B->coeffs, Blen, n);
    }

    _fmpz_poly_set_length(Q, n);
    _fmpz_poly_normalise(Q);
}
