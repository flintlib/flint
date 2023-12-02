/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void
_fmpq_poly_div_series(fmpz * Q, fmpz_t Qden,
                      const fmpz * A, const fmpz_t Aden, slong Alen,
                      const fmpz * B, const fmpz_t Bden, slong Blen, slong n)
{
    fmpz * C;
    fmpz_t Cden;

    C = _fmpz_vec_init(n);
    fmpz_init(Cden);

    _fmpq_poly_inv_series(C, Cden, B, Bden, Blen, n);
    _fmpq_poly_mullow(Q, Qden, A, Aden, Alen, C, Cden, n, n);

    _fmpz_vec_clear(C, n);
    fmpz_clear(Cden);
}

void fmpq_poly_div_series(fmpq_poly_t Q, const fmpq_poly_t A,
                                         const fmpq_poly_t B, slong n)
{
    if (A->length == 0)
    {
        fmpq_poly_zero(Q);
        return;
    }

    if (B->length == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_div_series). Division by zero.\n");
    }

    if (Q == A || Q == B)
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        fmpq_poly_div_series(t, A, B, n);
        fmpq_poly_swap(Q, t);
        fmpq_poly_clear(t);
        return;
    }

    fmpq_poly_fit_length(Q, n);
    _fmpq_poly_div_series(Q->coeffs, Q->den, A->coeffs, A->den, A->length,
                            B->coeffs, B->den, B->length, n);

    _fmpq_poly_set_length(Q, n);
    fmpq_poly_canonicalise(Q);
}

