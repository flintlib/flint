/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "fmpz.h"

void _fmpz_mod_poly_div_newton(fmpz * Q, const fmpz * A, slong lenA, 
                                     const fmpz * B, slong lenB, const fmpz_t p)
{
    const slong lenQ = lenA - lenB + 1;
    const slong Arev_len = lenQ + FLINT_MIN(lenB, lenQ);
    fmpz * Arev, * Brev;

    Arev = _fmpz_vec_init(Arev_len);
    Brev = Arev + lenQ;

    _fmpz_mod_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ);

    if (lenB >= lenQ)
        _fmpz_mod_poly_reverse(Brev, B + (lenB - lenQ), lenQ, lenQ);
    else
        _fmpz_mod_poly_reverse(Brev, B, lenB, lenB);

    _fmpz_mod_poly_div_series(Q, Arev, lenQ, Brev, lenB, p, lenQ);
    _fmpz_mod_poly_reverse(Q, Q, lenQ, lenQ);

    _fmpz_vec_clear(Arev, Arev_len);
}

void fmpz_mod_poly_div_newton(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                              const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length;
    const slong lenB = B->length;
    const slong lenQ = lenA - lenB + 1;

    fmpz * q;

    if (lenB == 0)
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
        {
            fmpz_mod_poly_set(Q, A, ctx);
            return;
        } else
        {
            flint_printf("Exception (fmpz_mod_poly_div_newton). Division by zero.\n");
            flint_abort();
        }
    }

    if (lenA < lenB)
    {
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    _fmpz_mod_poly_div_newton(q, A->coeffs, lenA, B->coeffs, lenB, fmpz_mod_ctx_modulus(ctx));

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenQ;
        Q->length = lenQ;
    } else
        _fmpz_mod_poly_set_length(Q, lenQ);
}
