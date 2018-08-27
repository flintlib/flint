/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nmod_mpoly.h"


void nmod_mpoly_pow_fmpz(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                  const fmpz_t pow, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * max_fields2;
    mp_bitcnt_t exp_bits;
    TMP_INIT;

    if (fmpz_cmp_ui(pow, UWORD(0)) < 0)
        flint_throw(FLINT_ERROR, "Negative power in nmod_mpoly_pow_fmpz");

    if (fmpz_fits_si(pow))
    {
        nmod_mpoly_pow_si(poly1, poly2, fmpz_get_si(pow), ctx);
        return;
    }

    /*
        we are raising a polynomial to an unreasonable exponent
        It must either be zero or a monomial with unit coefficient
    */

    if (poly2->length == WORD(0))
    {
        nmod_mpoly_zero(poly1, ctx);
        return;
    }

    if (poly2->length != WORD(1))
        flint_throw(FLINT_ERROR, "Multinomial in nmod_mpoly_pow_fmpz");

    TMP_START;

    max_fields2 = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(max_fields2 + i);
    mpoly_max_fields_fmpz(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);

    _fmpz_vec_scalar_mul_fmpz(max_fields2, max_fields2, ctx->minfo->nfields, pow);

    exp_bits = _fmpz_vec_max_bits(max_fields2, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    nmod_mpoly_fit_length(poly1, 1, ctx);
    nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
    poly1->bits = exp_bits;
    
    poly1->coeffs[0] = nmod_pow_fmpz(poly2->coeffs[0], pow, ctx->ffinfo->mod);
    mpoly_pack_vec_fmpz(poly1->exps + 0, max_fields2, exp_bits, ctx->minfo->nfields, 1);

    _nmod_mpoly_set_length(poly1, 1, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_clear(max_fields2 + i);

    TMP_END;
}
