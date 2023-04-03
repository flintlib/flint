/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr_poly.h"

void _fmpz_mod_poly_resultant_hgcd(fmpz_t res, const fmpz *A, slong lenA,
                               const fmpz *B, slong lenB, const fmpz_t mod)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz_mod(ctx, mod);  /* todo: by ref */
    GR_MUST_SUCCEED(_gr_poly_resultant_hgcd(res, A, lenA, B, lenB, FMPZ_MOD_POLY_HGCD_CUTOFF, FMPZ_MOD_POLY_GCD_CUTOFF, ctx));
    gr_ctx_clear(ctx);
}

void fmpz_mod_poly_resultant_hgcd(fmpz_t res, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    if (A->length == 0 || B->length == 0)
    {
       fmpz_zero(res);
    }
    else if (A->length < B->length)
    {
        fmpz_mod_poly_resultant_hgcd(res, B, A, ctx);

        if (((A->length | B->length) & 1) == 0)
           fmpz_negmod(res, res, fmpz_mod_ctx_modulus(ctx));
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length;

        /* lenA >= lenB >= 1 */
        _fmpz_mod_poly_resultant_hgcd(res, A->coeffs, lenA,
                                   B->coeffs, lenB, fmpz_mod_ctx_modulus(ctx));
    }
}
