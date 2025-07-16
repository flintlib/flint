/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"


/* if x = a/b, computes the cofactor of (b*x-a) in A */
void
_fmpz_poly_divexact_root_fmpq(fmpz * Q, const fmpz * A, slong len, const fmpq_t x)
{
    if (len < 2)
        return;

    fmpz_t r, t;

    fmpz_init(r);
    fmpz_init(t);

    // r = A[len - 1] / b
    fmpz_divexact(r, A + len - 1, fmpq_denref(x));

    for (slong i = len - 2; i > 0; i--)
    {
        // Q[i] = (A[i] + Q[i+1] * a) / b
        fmpz_set(t, A + i);
        fmpz_addmul(t, r, fmpq_numref(x));
        fmpz_divexact(t, t, fmpq_denref(x));

        fmpz_swap(Q + i, r);
        fmpz_swap(r, t);
    }
    // Q[0] = (A[0] + Q[1] * a) / b
    fmpz_swap(Q, r);

    fmpz_clear(r);
    fmpz_clear(t);
}

void
fmpz_poly_divexact_root_fmpq(fmpz_poly_t Q, const fmpz_poly_t A, const fmpq_t c)
{
    slong len = A->length;

    if (len <= 1)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (fmpq_is_zero(c))
    {
        fmpz_poly_shift_right(Q, A, 1);
        return;
    }

    fmpz_poly_fit_length(Q, len - 1);
    _fmpz_poly_divexact_root_fmpq(Q->coeffs, A->coeffs, len, c);
    _fmpz_poly_set_length(Q, len - 1);
}
