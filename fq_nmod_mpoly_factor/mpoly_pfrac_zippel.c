/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


slong fq_nmod_mpoly_set_eval_helper_and_zip_form2(
    slong * deg1_, /* degree of B wrt main var 1 */
    n_polyun_t EH,
    n_polyun_t H,
    n_polyun_t M,
    const fq_nmod_mpoly_t B,
    const fq_nmod_struct * betas,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong start, Bi, j, n;
    slong e0, e1, Hi, EHi;
    n_polyun_term_struct * EHterms, * Hterms, * Mterms;
    mp_limb_t * p;
    slong zip_length = 0;
    flint_bitcnt_t Bbits = B->bits;
    const fq_nmod_struct * Bcoeffs = B->coeffs;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Bbits);
    slong N = mpoly_words_per_exp_sp(Bbits, ctx->minfo);
    slong off0, off1, shift0, shift1;
    slong deg0, deg1 = -1;

    FLINT_ASSERT(Blen > 0);

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, Bbits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, Bbits, ctx->minfo);

    Bi = 0;
    deg0 = (Bexps[N*Bi + off0] >> shift0) & mask;

    EHi = 0;
    Hi = 0;

    while (Bi < Blen)
    {
        start = Bi;
        e0 = (Bexps[N*Bi + off0] >> shift0) & mask;
        e1 = (Bexps[N*Bi + off1] >> shift1) & mask;
        deg1 = FLINT_MAX(deg1, e1);
        while (1)
        {
            Bi++;
            if (Bi >= Blen)
                break;
            if (((Bexps[N*Bi + off0] >> shift0) & mask) != e0)
                break;
            if (((Bexps[N*Bi + off1] >> shift1) & mask) != e1)
                break;
        }

        n = Bi - start;

        n_polyun_fit_length(EH, EHi + 1);
        EHterms = EH->terms;
        EHterms[EHi].exp = pack_exp2(e0, e1);
        n_poly_fit_length(EHterms[EHi].coeff, d*3*n);
        EHterms[EHi].coeff->length = n;
        p = EHterms[EHi].coeff->coeffs;
        EHi++;

        _fq_nmod_mpoly_monomial_evals(p, Bexps + N*start, Bbits, n, betas, 2, ctx);

        if (e0 < deg0)
        {
            n_polyun_fit_length(H, Hi + 1);
            n_polyun_fit_length(M, Hi + 1);
            Hterms = H->terms;
            Mterms = M->terms;
            Hterms[Hi].exp = pack_exp2(e0, e1);
            Mterms[Hi].exp = pack_exp2(e0, e1);
            n_poly_fit_length(Hterms[Hi].coeff, d*n);
            zip_length = FLINT_MAX(zip_length, n);
            Hterms[Hi].coeff->length = n;
            flint_mpn_copyi(Hterms[Hi].coeff->coeffs, p, d*n);
            n_poly_fq_product_roots_n_fq(Mterms[Hi].coeff, p, n, ctx->fqctx);
            Hi++;
        }

        for (j = n - 1; j >= 0; j--)
        {
            _n_fq_set(p + d*(3*j + 2), p + d*j, d);
            _n_fq_set(p + d*(3*j + 0), p + d*(3*j + 2), d);
            n_fq_set_fq_nmod(p + d*(3*j + 1), Bcoeffs + start + j, ctx->fqctx);
        }
    }

    EH->length = EHi;
    H->length = Hi;
    M->length = Hi;

    *deg1_ = deg1;
    return zip_length;
}
