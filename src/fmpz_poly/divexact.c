/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "gr.h"
#include "gr_poly.h"

void
_fmpz_poly_divexact(fmpz * Q, const fmpz * A, slong lenA,
                                         const fmpz * B, slong lenB)
{
    slong lenQ = lenA - lenB + 1;

    if (lenQ == 1)
    {
        fmpz_divexact(Q, A + lenA - 1, B + lenB - 1);
    }
    else if (lenB == 1)
    {
        if (fmpz_is_pm1(B))
            _fmpz_vec_scalar_mul_fmpz(Q, A, lenA, B);
        else
            _fmpz_vec_scalar_divexact_fmpz(Q, A, lenA, B);
    }
    else if (lenQ <= 100 || lenB <= 16)
    {
        gr_ctx_t ctx;
        gr_ctx_init_fmpz(ctx);
        GR_MUST_SUCCEED(_gr_poly_divexact_basecase_bidirectional(Q, A, lenA, B, lenB, ctx));
    }
    else
    {
        /* todo: the true cutoffs are sometimes higher, especially with
           unbalanced operands, possibly because divconquer division needs tuning */
        slong A_bits, B_bits, B_cutoff, Q_cutoff;
        gr_ctx_t ctx;
        gr_ctx_init_fmpz(ctx);

        A_bits = _fmpz_vec_max_bits(A, lenQ);
        B_bits = _fmpz_vec_max_bits(B, FLINT_MIN(lenB, lenQ));
        A_bits = FLINT_ABS(A_bits);
        B_bits = FLINT_ABS(B_bits);

        B_cutoff = (B_bits > 3000) ? 20 : 60;
        Q_cutoff = (A_bits > 1000) ? 100 : 200;

        if (A_bits >= 100 * B_bits)
        {
            Q_cutoff *= 2;
            B_cutoff *= 2;
        }

        if (lenQ <= Q_cutoff || lenB <= B_cutoff)
            GR_MUST_SUCCEED(_gr_poly_divexact_basecase_bidirectional(Q, A, lenA, B, lenB, ctx));
        else
            _fmpz_poly_div(Q, A, lenA, B, lenB, 0);
    }

}

void
fmpz_poly_divexact(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
    fmpz_poly_t T;
    slong lenA = A->length;
    slong lenB = B->length;
    slong lenQ = lenA - lenB + 1;

    if (lenB == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_divexact). Division by zero.\n");
    }

    if (lenA < lenB)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpz_poly_init2(T, lenQ);
        _fmpz_poly_divexact(T->coeffs, A->coeffs, lenA, B->coeffs, lenB);
        _fmpz_poly_set_length(T, lenQ);
        fmpz_poly_swap(T, Q);
        fmpz_poly_clear(T);
    }
    else
    {
        fmpz_poly_fit_length(Q, lenQ);
        _fmpz_poly_divexact(Q->coeffs, A->coeffs, lenA, B->coeffs, lenB);
        _fmpz_poly_set_length(Q, lenQ);
    }

    /* should not be needed, but produce something normalised in case
       this was called with invalid input */
    _fmpz_poly_normalise(Q);
}
