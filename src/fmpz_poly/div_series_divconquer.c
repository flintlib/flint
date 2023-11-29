/*
    Copyright (C) 2019 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_div_series_divconquer(fmpz * Q, const fmpz * A, slong Alen,
    const fmpz * B, slong Blen, slong n)
{
    fmpz * Arev = (fmpz *) _fmpz_vec_init(2*n - 1);
    fmpz * Brev = (fmpz *) _fmpz_vec_init(n);

    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);

    _fmpz_poly_reverse(Arev, A, Alen, 2*n - 1);
    _fmpz_poly_reverse(Brev, B, Blen, n);

    if (!_fmpz_poly_div(Q, Arev, 2*n - 1, Brev, n, 1))
    {
        _fmpz_vec_clear(Arev, 2*n - 1); /* flint_throw */
        _fmpz_vec_clear(Brev, n); /* flint_throw */
        flint_throw(FLINT_ERROR, "Not an exact division\n");
    }

    _fmpz_poly_reverse(Q, Q, n, n);

    _fmpz_vec_clear(Arev, 2*n - 1);
    _fmpz_vec_clear(Brev, n);
}

void fmpz_poly_div_series_divconquer(fmpz_poly_t Q, const fmpz_poly_t A,
                                         const fmpz_poly_t B, slong n)
{
    slong Alen = FLINT_MIN(A->length, n);
    slong Blen = FLINT_MIN(B->length, n);

    if (Blen == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_div_series_divconquer). Division by zero.\n");
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
        _fmpz_poly_div_series_divconquer(t->coeffs, A->coeffs, Alen, B->coeffs, Blen, n);
        fmpz_poly_swap(Q, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(Q, n);
        _fmpz_poly_div_series_divconquer(Q->coeffs, A->coeffs, Alen, B->coeffs, Blen, n);
    }

    _fmpz_poly_set_length(Q, n);
    _fmpz_poly_normalise(Q);
}

