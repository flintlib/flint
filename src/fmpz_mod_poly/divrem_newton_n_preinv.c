/*
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_divrem_newton_n_preinv (fmpz* Q, fmpz* R, const fmpz* A,
                              slong lenA, const fmpz* B, slong lenB,
                              const fmpz* Binv, slong lenBinv, const fmpz_mod_ctx_t ctx)
{
    const slong lenQ = lenA - lenB + 1;

    _fmpz_mod_poly_div_newton_n_preinv(Q, A, lenA, B, lenB, Binv, lenBinv, ctx);

    if (lenB > 1)
    {
        if (lenQ >= lenB - 1)
            _fmpz_mod_poly_mullow(R, Q, lenQ, B, lenB - 1, lenB - 1, ctx);
        else
            _fmpz_mod_poly_mullow(R, B, lenB - 1, Q, lenQ, lenB - 1, ctx);

        _fmpz_mod_vec_sub(R, A, R, lenB - 1, ctx);
    }
}

void fmpz_mod_poly_divrem_newton_n_preinv(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
                              const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                         const fmpz_mod_poly_t Binv, const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length, lenB = B->length, lenBinv= Binv->length,
                       lenQ = lenA - lenB + 1;
    fmpz* q, * r;

    if (lenB == 0)
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
        {
            fmpz_mod_poly_set(Q, A, ctx);
            fmpz_mod_poly_zero(R, ctx);
            return;
        }
        else
        {
            flint_throw(FLINT_ERROR, "(fmpz_mod_poly_divrem_newton_n_preinv): Division by zero.\n");
        }
    }

    if (lenA < lenB)
    {
        fmpz_mod_poly_set(R, A, ctx);
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

    if (lenA > 2 * lenB - 2)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_divrem_newton_n_preinv).\n");
    }

    if (Q == A || Q == B || Q == Binv)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }
    if (R == A || R == B || R == Binv)
    {
        r = _fmpz_vec_init(lenB - 1);
    }
    else
    {
        fmpz_mod_poly_fit_length(R, lenB - 1, ctx);
        r = R->coeffs;
    }

    _fmpz_mod_poly_divrem_newton_n_preinv (q, r, A->coeffs, lenA,
                                           B->coeffs, lenB, Binv->coeffs,
                                           lenBinv, ctx);

    if (Q == A || Q == B || Q == Binv)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _fmpz_mod_poly_set_length(Q, lenQ);
    }
    if (R == A || R == B || R == Binv)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenB - 1;
        R->length = lenB - 1;
    }
    else
    {
      _fmpz_mod_poly_set_length(R, lenB - 1);
    }

    _fmpz_mod_poly_normalise(R);
}
