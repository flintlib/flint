/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "n_poly.h"

/*
    only E and alphas are shifted by "var"
    so output is in E[0]
    and first relevent alpha is alphas[0]
*/
int _fmpz_mpoly_evaluate_rest_fmpz(
    fmpz * E,
    slong * starts,
    slong * ends,
    slong * stops,
    ulong * es,
    const fmpz * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    slong var,
    const fmpz * alphas,
    const slong * offsets,
    const slong * shifts,
    slong N,
    ulong mask,
    slong nvars)
{
    slong v, stop;
    ulong next_e;

    FLINT_ASSERT(var < nvars);

    E -= var;
    alphas -= var;

    v = var;
    starts[v] = 0;
    ends[v] = Alen;
    fmpz_zero(E + v);

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

    fmpz_zero(E + v);

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
        fmpz_add(E + v, E + v, E + v + 1);
    }
    else
    {
        fmpz_add(E + v, E + v, Acoeffs + starts[v]);
    }

    if (stops[v] < ends[v])
    {
        next_e = mask & (Aexps[N*stops[v] + offsets[v]] >> shifts[v]);
        FLINT_ASSERT(next_e < es[v]);
        fmpz_pow_ui(E + v + 1, alphas + v, es[v] - next_e);
        fmpz_mul(E + v, E + v, E + v + 1);
        es[v] = next_e;
        starts[v] = stops[v];
        goto next;
    }
    else
    {
        fmpz_pow_ui(E + v + 1, alphas + v, es[v]);
        fmpz_mul(E + v, E + v, E + v + 1);
    }

    if (v > var)
    {
        v--;
        goto calculate_return;
    }

    return 1;
}


void _fmpz_mpoly_eval_rest_to_poly(
    fmpz_poly_t E,
    const fmpz_mpoly_t A,
    const fmpz * alphas,
    const fmpz_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong * offsets, * shifts;
    slong offset, shift;
    slong start, stop;
    ulong e, mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong * starts, * ends, * stops;
    ulong * es;
    fmpz * realE;

    FLINT_ASSERT(n > 1);

    E->length = 0;
    if (A->length < 1)
        return;

    starts = FLINT_ARRAY_ALLOC(n, slong);
    ends   = FLINT_ARRAY_ALLOC(n, slong);
    stops  = FLINT_ARRAY_ALLOC(n, slong);
    es     = FLINT_ARRAY_ALLOC(n, ulong);
    realE  = FLINT_ARRAY_ALLOC(n + 1, fmpz);
    for (i = 0; i < n + 1; i++)
        fmpz_init(realE + i);

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

    fmpz_poly_fit_length(E, e + 1);
    while (E->length <= e)
    {
        fmpz_zero(E->coeffs + E->length);
        E->length++;
    }

    _fmpz_mpoly_evaluate_rest_fmpz(realE, starts, ends, stops, es,
                    A->coeffs + start, A->exps + N*start, stop - start, 1,
                          alphas, offsets, shifts, N, mask, ctx->minfo->nvars);
    fmpz_set(E->coeffs + e, realE + 0);

    if (stop < A->length)
    {
        FLINT_ASSERT(e > (mask & (A->exps[N*stop + offset] >> shift)));
        e = (mask & (A->exps[N*stop + offset] >> shift));
        start = stop;
        goto next;
    }

    _fmpz_poly_normalise(E);

    for (i = 0; i < n + 1; i++)
        fmpz_clear(realE + i);
    flint_free(realE);
    flint_free(starts);
    flint_free(ends);
    flint_free(stops);
    flint_free(es);

    flint_free(offsets);
    flint_free(shifts);
}

