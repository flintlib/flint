/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


void nmod_mpoly_get_coeff_vars_ui(nmod_mpoly_t C, const nmod_mpoly_t A,
                         const slong * vars, const ulong * exps, slong length,
                                                    const nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong i, j;
    slong offset, shift;
    slong maxoffset, minoffset;
    ulong * uexp;
    ulong * tmask, * texp;
    slong nvars = ctx->minfo->nvars;
    mp_limb_t * Ccoeff;
    ulong * Cexp;
    slong Clen;
    TMP_INIT;

    if (C == A)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init(T, ctx);
        nmod_mpoly_get_coeff_vars_ui(T, A, vars, exps,length, ctx);
        nmod_mpoly_swap(T, C, ctx);
        nmod_mpoly_clear(T, ctx);
        return;
    }

    TMP_START;

    uexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    for (i = 0; i < nvars; i++)
        uexp[i] = 0;
    for (i = 0; i < length; i++)
        uexp[vars[i]] = exps[i];

    if (bits < mpoly_exp_bits_required_ui(uexp, ctx->minfo))
    {
        nmod_mpoly_zero(C, ctx);
        goto cleanup;
    }

    nmod_mpoly_fit_length_reset_bits(C, 4, bits, ctx);

    tmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(tmask, N);
    mpoly_set_monomial_ui(texp, uexp, bits, ctx->minfo);

    if (bits <= FLINT_BITS)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        maxoffset = 0;
        minoffset = N;
        for (i = 0; i < length; i++)
        {
            mpoly_gen_offset_shift_sp(&offset, &shift, vars[i], bits, ctx->minfo);
            tmask[offset] |= mask << shift;
            maxoffset = FLINT_MAX(maxoffset, offset);
            minoffset = FLINT_MIN(minoffset, offset);
        }
        FLINT_ASSERT(minoffset < N);

        Ccoeff = C->coeffs;
        Cexp = C->exps;
        Clen = 0;
        for (i = 0; i < A->length; i++)
        {
            for (j = minoffset; j <= maxoffset; j++)
            {
                if ((((A->exps + N*i)[j] ^ texp[j]) & tmask[j]) != UWORD(0))
                    goto continue_outer_sp;
            }
            _nmod_mpoly_fit_length(&Ccoeff, &C->coeffs_alloc,
                                   &Cexp, &C->exps_alloc, N, Clen + 1);
            mpoly_monomial_sub(Cexp + N*Clen, A->exps + N*i, texp, N);
            Ccoeff[Clen] = A->coeffs[i];
            Clen++;
continue_outer_sp:;
        }

        C->coeffs = Ccoeff;
        C->exps = Cexp;
        _nmod_mpoly_set_length(C, Clen, ctx);
    }
    else
    {
        ulong wpf = A->bits/FLINT_BITS;
        maxoffset = 0;
        minoffset = N;
        for (i = 0; i < length; i++)
        {
            offset = mpoly_gen_offset_mp(vars[i], A->bits, ctx->minfo);
            minoffset = FLINT_MIN(minoffset, offset);
            maxoffset = FLINT_MAX(maxoffset, offset + wpf - 1);
            for (j = 0; j < wpf; j++)
                tmask[offset + j] = -UWORD(1);
        }
        FLINT_ASSERT(minoffset < N);

        Ccoeff = C->coeffs;
        Cexp = C->exps;
        Clen = 0;
        for (i = 0; i < A->length; i++)
        {
            for (j = minoffset; j <= maxoffset; j++)
            {
                if ((((A->exps + N*i)[j] ^ texp[j]) & tmask[j]) != UWORD(0))
                    goto continue_outer_mp;
            }
            _nmod_mpoly_fit_length(&Ccoeff, &C->coeffs_alloc,
                                   &Cexp, &C->exps_alloc, N, Clen + 1);
            mpoly_monomial_sub_mp(Cexp + N*Clen, A->exps + N*i, texp, N);
            Ccoeff[Clen] = A->coeffs[i];
            Clen++;
continue_outer_mp:;
        }

        C->coeffs = Ccoeff;
        C->exps = Cexp;
        _nmod_mpoly_set_length(C, Clen, ctx);
    }

cleanup:

    TMP_END;
    return;
}
