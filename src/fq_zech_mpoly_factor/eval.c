/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"


/*
    only E and alphas are shifted by "var"
    so output is in E[0]
*/
int _fq_zech_mpoly_eval_rest_fq_zech_poly(
    fq_zech_poly_struct * E,
    slong * starts,
    slong * ends,
    slong * stops,
    ulong * es,
    const fq_zech_struct * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    slong var,
    const fq_zech_poly_struct * alphas,
    const slong * offsets,
    const slong * shifts,
    slong N,
    ulong mask,
    slong nvars,
    const fq_zech_ctx_t ctx)
{
    slong v, stop;
    ulong next_e;

    FLINT_ASSERT(var < nvars);

    E -= var;
    alphas -= var;

    v = var;
    starts[v] = 0;
    ends[v] = Alen;
    fq_zech_poly_zero(E + v, ctx);

    if (Alen < 1)
        return 1;

calculate:
/*
    input:
        v
        starts[v]
        ends[v]
*/
    FLINT_ASSERT(ends[v] > starts[v]);
    es[v] = mask & (Aexps[N*starts[v] + offsets[v]] >> shifts[v]);

    fq_zech_poly_zero(E + v, ctx);

next:

    FLINT_ASSERT(es[v] == (mask & (Aexps[N*starts[v] + offsets[v]] >> shifts[v])));

    stop = starts[v] + 1;
    while (stop < ends[v] &&
           (mask & (Aexps[N*stop + offsets[v]] >> shifts[v])) == es[v])
    {
        stop++;
    }
    stops[v] = stop;

    if (v + 1 < nvars)
    {
        starts[v + 1] = starts[v];
        ends[v + 1] = stops[v];
        v++;
        goto calculate;
calculate_return:
        fq_zech_poly_add(E + v, E + v, E + v + 1, ctx);
    }
    else
    {
        fq_zech_poly_set_fq_zech(E + v + 1, Acoeffs + starts[v], ctx);
        fq_zech_poly_add(E + v, E + v, E + v + 1, ctx);
    }

    if (stops[v] < ends[v])
    {
        next_e = mask & (Aexps[N*stops[v] + offsets[v]] >> shifts[v]);
        FLINT_ASSERT(next_e < es[v]);
        fq_zech_poly_pow(E + v + 1, alphas + v, es[v] - next_e, ctx);
        fq_zech_poly_mul(E + v, E + v, E + v + 1, ctx);
        es[v] = next_e;
        starts[v] = stops[v];
        goto next;
    }
    else
    {
        fq_zech_poly_pow(E + v + 1, alphas + v, es[v], ctx);
        fq_zech_poly_mul(E + v, E + v, E + v + 1, ctx);
    }

    if (v > var)
    {
        v--;
        goto calculate_return;
    }

    return 1;
}


void _fq_zech_mpoly_eval_to_bpoly(
    fq_zech_bpoly_t E,
    const fq_zech_mpoly_t A,
    const fq_zech_poly_struct * alphabetas,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong * offsets, * shifts;
    slong offset, shift;
    slong start, stop;
    ulong e, mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong * starts, * ends, * stops;
    ulong * es;
    fq_zech_poly_struct * realE;

    E->length = 0;
    if (A->length < 1)
        return;

    starts = FLINT_ARRAY_ALLOC(n, slong);
    ends   = FLINT_ARRAY_ALLOC(n, slong);
    stops  = FLINT_ARRAY_ALLOC(n, slong);
    es     = FLINT_ARRAY_ALLOC(n, ulong);
    realE  = FLINT_ARRAY_ALLOC(n + 1, fq_zech_poly_struct);
    for (i = 0; i < n + 1; i++)
        fq_zech_poly_init(realE + i, ctx->fqctx);

    offsets = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);
    shifts  = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);
    for (i = 0; i < ctx->minfo->nvars; i++)
        mpoly_gen_offset_shift_sp(offsets + i, shifts + i, i, A->bits, ctx->minfo);

    offset = offsets[0];
    shift = shifts[0];

    start = 0;
    e = mask & (A->exps[N*start + offset] >> shift);

next:

    FLINT_ASSERT(start < A->length);
    FLINT_ASSERT(e == (mask & (A->exps[N*start + offset] >> shift)));

    stop = start + 1;
    while (stop < A->length && (mask & (A->exps[N*stop + offset] >> shift)) == e)
        stop++;

    fq_zech_bpoly_fit_length(E, e + 1, ctx->fqctx);
    while (E->length <= e)
    {
        fq_zech_poly_zero(E->coeffs + E->length, ctx->fqctx);
        E->length++;
    }

    _fq_zech_mpoly_eval_rest_fq_zech_poly(realE, starts, ends, stops, es,
                    A->coeffs + start, A->exps + N*start, stop - start, 1,
                                        alphabetas, offsets, shifts, N, mask,
                                                ctx->minfo->nvars, ctx->fqctx);
    fq_zech_poly_set(E->coeffs + e, realE + 0, ctx->fqctx);

    if (stop < A->length)
    {
        FLINT_ASSERT(e > (mask & (A->exps[N*stop + offset] >> shift)));
        e = (mask & (A->exps[N*stop + offset] >> shift));
        start = stop;
        goto next;
    }

    fq_zech_bpoly_normalise(E, ctx->fqctx);

    for (i = 0; i < n + 1; i++)
        fq_zech_poly_clear(realE + i, ctx->fqctx);
    flint_free(realE);
    flint_free(starts);
    flint_free(ends);
    flint_free(stops);
    flint_free(es);

    flint_free(offsets);
    flint_free(shifts);
}


/* A = B(gen(var), 0) */
void _fq_zech_mpoly_set_fq_zech_bpoly_var1_zero(
    fq_zech_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_zech_bpoly_t B,
    slong var,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);
    slong i, Alen;
    slong Blen = B->length;
    ulong * genexp;
    TMP_INIT;

    TMP_START;

    genexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    if (Abits <= FLINT_BITS)
        mpoly_gen_monomial_sp(genexp, var, Abits, ctx->minfo);
    else
        mpoly_gen_monomial_offset_mp(genexp, var, Abits, ctx->minfo);

    Alen = 2;
    for (i = 0; i < Blen; i++)
        Alen += (B->coeffs[i].length > 0);

    fq_zech_mpoly_fit_length_reset_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        FLINT_ASSERT(Alen < A->alloc);
        fq_zech_bpoly_get_coeff(A->coeffs + Alen, B, i, 0, ctx->fqctx);
        if (fq_zech_is_zero(A->coeffs + Alen, ctx->fqctx))
            continue;

        if (Abits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps + N*Alen, genexp, N, i);
        else
            mpoly_monomial_mul_ui_mp(A->exps + N*Alen, genexp, N, i);
        Alen++;
    }
    A->length = Alen;

    TMP_END;
}
