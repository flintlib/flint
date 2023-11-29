/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_div_newton_n_preinv (fmpz* Q, const fmpz* A, slong lenA,
                                   const fmpz* B, slong lenB, const fmpz* Binv,
                                   slong lenBinv, const fmpz_mod_ctx_t ctx)
{
    const slong lenQ = lenA - lenB + 1;
    fmpz * Arev;

    Arev = _fmpz_vec_init(lenQ);

    _fmpz_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ);
    _fmpz_mod_poly_mullow(Q, Arev, lenQ, Binv, FLINT_MIN(lenQ, lenBinv), lenQ, ctx);
    _fmpz_poly_reverse(Q, Q, lenQ, lenQ);

    _fmpz_vec_clear(Arev, lenQ);
}


void fmpz_mod_poly_div_newton_n_preinv(fmpz_mod_poly_t Q,
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                          const fmpz_mod_poly_t Binv, const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1,
                lenBinv = Binv->length;
    fmpz *q;

    if (lenB == 0)
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
        {
            fmpz_mod_poly_set(Q, A, ctx);
            return;
        }
        else
        {
            flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_div_newton_n_preinv). Division by zero.\n");
        }
    }

    if (lenA < lenB)
    {
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

    if (lenA > 2 * lenB - 2)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_div_newton_n_preinv).\n");
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

    _fmpz_mod_poly_div_newton_n_preinv (q, A->coeffs, lenA, B->coeffs, lenB,
                             Binv->coeffs, lenBinv, ctx);

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
}
