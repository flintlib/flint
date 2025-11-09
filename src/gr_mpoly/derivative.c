/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Based on src/fmpz_mpoly/derivative.c */

#include "fmpz.h"
#include "mpoly.h"
#include "gr_mpoly.h"

static int
_gr_mpoly_derivative_sp(gr_ptr coeff1, ulong * exp1, slong * len1,
                        gr_srcptr coeff2, const ulong * exp2, slong len2,
                        flint_bitcnt_t bits, slong N, slong offset, slong shift,
                        ulong * oneexp, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;
    slong i;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);

    /* x^c -> c*x^(c-1) */
    *len1 = 0;
    for (i = 0; i < len2; i++)
    {
        ulong c = (exp2[N*i + offset] >> shift) & mask;
        gr_ptr c1 = GR_ENTRY(coeff1, *len1, cctx->sizeof_elem);
        status |= gr_mul_ui(c1, GR_ENTRY(coeff2, i, cctx->sizeof_elem), c, cctx);
        if (gr_is_zero(c1, cctx) != T_TRUE)
        {
            mpoly_monomial_sub(exp1 + N*(*len1), exp2 + N*i, oneexp, N);
            (*len1)++;
        }
    }

    return status;
}


static int
_gr_mpoly_derivative_mp(gr_ptr coeff1, ulong * exp1, slong * len1,
                        gr_srcptr coeff2, const ulong * exp2, slong len2,
                        flint_bitcnt_t bits, slong N, slong offset,
                        ulong * oneexp, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;
    slong i;
    fmpz_t c;
    fmpz_init(c);

    /* x^c -> c*x^(c-1) */
    *len1 = 0;
    for (i = 0; i < len2; i++)
    {
        fmpz_set_ui_array(c, exp2 + N*i + offset, bits/FLINT_BITS);
        gr_ptr c1 = GR_ENTRY(coeff1, *len1, cctx->sizeof_elem);
        status |= gr_mul_fmpz(c1, GR_ENTRY(coeff2, i, cctx->sizeof_elem), c, cctx);
        if (gr_is_zero(c1, cctx) != T_TRUE)
        {
            mpoly_monomial_sub_mp(exp1 + N*(*len1), exp2 + N*i, oneexp, N);
            (*len1)++;
        }
    }

    fmpz_clear(c);

    return status;
}


int
gr_mpoly_derivative(gr_mpoly_t A, const gr_mpoly_t B, slong var,
                    gr_mpoly_ctx_t ctx)
{
    int status = GR_SUCCESS;
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong N, offset, shift;
    flint_bitcnt_t bits;
    ulong * oneexp;
    slong len1;
    TMP_INIT;

    TMP_START;
    bits = B->bits;

    if (A != B)
        gr_mpoly_fit_length_reset_bits(A, B->length, bits, ctx);

    N = mpoly_words_per_exp(bits, GR_MPOLY_MCTX(ctx));
    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    if (bits <= FLINT_BITS)
    {
        mpoly_gen_monomial_offset_shift_sp(oneexp, &offset, &shift,
                                           var, bits, mctx);
        status |= _gr_mpoly_derivative_sp(A->coeffs, A->exps, &len1,
                                          B->coeffs, B->exps, B->length,
                                          bits, N, offset, shift, oneexp, cctx);
    }
    else
    {
        offset = mpoly_gen_monomial_offset_mp(oneexp, var, bits, mctx);
        status |= _gr_mpoly_derivative_mp(A->coeffs, A->exps, &len1,
                                          B->coeffs, B->exps, B->length,
                                          bits, N, offset, oneexp, cctx);
    }

    _gr_mpoly_set_length(A, len1, ctx);

    TMP_END;

    return status;
}
