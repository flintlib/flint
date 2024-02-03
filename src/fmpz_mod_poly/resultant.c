/*
    Copyright (C) 2011, 2014 William Hart
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr_poly.h"

void _fmpz_mod_poly_resultant(fmpz_t res, const fmpz *A, slong lenA,
                               const fmpz *B, slong lenB, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_resultant(res, A, lenA, B, lenB, gr_ctx));
}

void fmpz_mod_poly_resultant(fmpz_t res, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    if (A->length == 0 || B->length == 0)
    {
        fmpz_zero(res);
    }
    else if (A->length < B->length)
    {
        fmpz_mod_poly_resultant(res, B, A, ctx);

        if (((A->length | B->length) & 1) == 0)
            fmpz_mod_neg(res, res, ctx);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length;

        /* lenA >= lenB >= 1 */
        _fmpz_mod_poly_resultant(res, A->coeffs, lenA, B->coeffs, lenB, ctx);
    }
}
