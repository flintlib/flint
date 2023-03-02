/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

/*
    only E and alphas are shifted by "var"
    so output is in E[0]
*/
int _nmod_mpoly_evaluate_rest_n_poly(
    n_poly_struct * E,
    slong * starts,
    slong * ends,
    slong * stops,
    ulong * es,
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    slong var,
    const n_poly_struct * alphas,
    const slong * offsets,
    const slong * shifts,
    slong N,
    ulong mask,
    slong nvars,
    nmod_t ctx)
{
    slong v, stop;
    ulong next_e;

    FLINT_ASSERT(var < nvars);

    E -= var;
    alphas -= var;

    v = var;
    starts[v] = 0;
    ends[v] = Alen;
    n_poly_zero(E + v);

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

    n_poly_zero(E + v);

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
        n_poly_mod_add(E + v, E + v, E + v + 1, ctx);
    }
    else
    {
        n_poly_mod_add_ui(E + v, E + v, Acoeffs[starts[v]], ctx);
    }

    if (stops[v] < ends[v])
    {
        next_e = mask & (Aexps[N*stops[v] + offsets[v]] >> shifts[v]);
        FLINT_ASSERT(next_e < es[v]);
        n_poly_mod_pow(E + v + 1, alphas + v, es[v] - next_e, ctx);
        n_poly_mod_mul(E + v, E + v, E + v + 1, ctx);
        es[v] = next_e;
        starts[v] = stops[v];
        goto next;
    }
    else
    {
        n_poly_mod_pow(E + v + 1, alphas + v, es[v], ctx);
        n_poly_mod_mul(E + v, E + v, E + v + 1, ctx);
    }

    if (v > var)
    {
        v--;
        goto calculate_return;
    }

    return 1;
}


void _nmod_mpoly_eval_rest_to_n_bpoly(
    n_bpoly_t E,
    const nmod_mpoly_t A,
    const n_poly_struct * alphabetas,
    const nmod_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong * offsets, * shifts;
    slong offset, shift;
    slong start, stop;
    ulong e, mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong * starts, * ends, * stops;
    ulong * es;
    n_poly_struct * realE;

    E->length = 0;
    if (A->length < 1)
        return;

    starts = FLINT_ARRAY_ALLOC(n, slong);
    ends   = FLINT_ARRAY_ALLOC(n, slong);
    stops  = FLINT_ARRAY_ALLOC(n, slong);
    es     = FLINT_ARRAY_ALLOC(n, ulong);
    realE  = FLINT_ARRAY_ALLOC(n + 1, n_poly_struct);
    for (i = 0; i < n + 1; i++)
        n_poly_init(realE + i);

    offsets = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
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

    n_bpoly_fit_length(E, e + 1);
    while (E->length <= e)
    {
        n_poly_zero(E->coeffs + E->length);
        E->length++;
    }

    _nmod_mpoly_evaluate_rest_n_poly(realE, starts, ends, stops, es,
                    A->coeffs + start, A->exps + N*start, stop - start, 1,
                                        alphabetas, offsets, shifts, N, mask,
                                          ctx->minfo->nvars, ctx->mod);
    n_poly_set(E->coeffs + e, realE + 0);

    if (stop < A->length)
    {
        FLINT_ASSERT(e > (mask & (A->exps[N*stop + offset] >> shift)));
        e = (mask & (A->exps[N*stop + offset] >> shift));
        start = stop;
        goto next;
    }

    n_bpoly_normalise(E);

    for (i = 0; i < n + 1; i++)
        n_poly_clear(realE + i);
    flint_free(realE);
    flint_free(starts);
    flint_free(ends);
    flint_free(stops);
    flint_free(es);

    flint_free(offsets);
    flint_free(shifts);
}

/* A = B(gen(var), 0) */
void _nmod_mpoly_set_n_bpoly_var1_zero(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const n_bpoly_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
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

    nmod_mpoly_fit_length_reset_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        mp_limb_t c = n_poly_get_coeff(B->coeffs + i, 0);
        if (c == 0)
            continue;

        FLINT_ASSERT(Alen < A->coeffs_alloc);
        A->coeffs[Alen] = c;
        if (Abits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps + N*Alen, genexp, N, i);
        else
            mpoly_monomial_mul_ui_mp(A->exps + N*Alen, genexp, N, i);
        Alen++;
    }
    A->length = Alen;

    TMP_END;
}

