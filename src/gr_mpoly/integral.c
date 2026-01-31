/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Based on src/fmpz_mpoly/integral.c */

#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpoly.h"
#include "gr_mpoly.h"

static int
_gr_mpoly_integral(gr_ptr coeff1, ulong * exp1,
                   gr_srcptr coeff2, const ulong * exp2, slong len2,
                   slong var, slong bits, gr_mpoly_ctx_t ctx)
{
    int status = GR_SUCCESS;
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong i, N;
    slong offset;
    slong shift;
    ulong * oneexp;
    TMP_INIT;

    TMP_START;

    N = mpoly_words_per_exp(bits, mctx);
    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    if (bits <= FLINT_BITS)
    {
        ulong c, mask;

        mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        mpoly_gen_monomial_offset_shift_sp(oneexp, &offset, &shift, var, bits, mctx);

        /* x^c -> x^(c+1)/(c+1) */
        for (i = 0; i < len2; i++)
        {
            c = (exp2[N*i + offset] >> shift) & mask;
            c++;
            gr_ptr c1 = GR_ENTRY(coeff1, i, cctx->sizeof_elem);
            status |= gr_div_ui(c1, GR_ENTRY(coeff2, i, cctx->sizeof_elem), c, cctx);
            mpoly_monomial_add(exp1 + N*i, exp2 + N*i, oneexp, N);
        }
    }
    else
    {
        fmpz_t c;
        fmpz_init(c);

        offset = mpoly_gen_monomial_offset_mp(oneexp, var, bits, mctx);

        /* x^c -> x^(c+1)/(c+1) */
        for (i = 0; i < len2; i++)
        {
            fmpz_set_ui_array(c, exp2 + N*i + offset, bits/FLINT_BITS);
            fmpz_add_ui(c, c, 1);
            gr_ptr c1 = GR_ENTRY(coeff1, i, cctx->sizeof_elem);
            status |= gr_div_fmpz(c1, GR_ENTRY(coeff2, i, cctx->sizeof_elem), c, cctx);
            mpoly_monomial_add_mp(exp1 + N*i, exp2 + N*i, oneexp, N);
        }

        fmpz_clear(c);
    }

    TMP_END;

    return status;
}


int
gr_mpoly_integral(gr_mpoly_t poly1, const gr_mpoly_t poly2, slong var,
                  gr_mpoly_ctx_t ctx)
{
    int status = GR_SUCCESS;
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    slong i;
    flint_bitcnt_t exp_bits;
    ulong * exp2 = poly2->exps;
    fmpz * gen_fields, * max_fields;
    int free2 = 0;
    TMP_INIT;

    TMP_START;

    /* compute bits required to represent result */
    gen_fields = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    max_fields = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
    {
        fmpz_init(gen_fields + i);
        fmpz_init(max_fields + i);
    }
    mpoly_gen_fields_fmpz(gen_fields, var, mctx);
    mpoly_max_fields_fmpz(max_fields, poly2->exps, poly2->length,
                          poly2->bits, mctx);
    _fmpz_vec_add(max_fields, max_fields, gen_fields, mctx->nfields);

    exp_bits = _fmpz_vec_max_bits(max_fields, mctx->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);

    for (i = 0; i < mctx->nfields; i++)
    {
        fmpz_clear(gen_fields + i);
        fmpz_clear(max_fields + i);
    }


    /* ensure input exponents are packed into same sized fields as output */
    if (exp_bits > poly2->bits)
    {
        slong N = mpoly_words_per_exp(exp_bits, mctx);
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                               poly2->length, mctx);
    }

    /* deal with aliasing and do integration */
    if (poly1 == poly2)
    {
        gr_mpoly_t temp;

        gr_mpoly_init2(temp, poly2->length, ctx);
        gr_mpoly_fit_bits(temp, exp_bits, ctx);
        temp->bits = exp_bits;

        status |= _gr_mpoly_integral(temp->coeffs, temp->exps,
                                     poly2->coeffs, exp2, poly2->length,
                                     var, exp_bits, ctx);
        _gr_mpoly_set_length(temp, poly2->length, ctx);

        gr_mpoly_swap(temp, poly1, ctx);
        gr_mpoly_clear(temp, ctx);
    }
    else
    {
        gr_mpoly_fit_length(poly1, poly2->length, ctx);
        gr_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        status |= _gr_mpoly_integral(poly1->coeffs, poly1->exps,
                                     poly2->coeffs, exp2, poly2->length,
                                     var, exp_bits, ctx);
        _gr_mpoly_set_length(poly1, poly2->length, ctx);
    }

    if (free2)
        flint_free(exp2);

    TMP_END;

    return status;
}
