/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"
#include "fq_nmod_mpoly_factor.h"

/*
    content wrt gen(0), ..., gen(num_vars-1)
    successful answer will be returned with g->bits == A->bits
*/
int fq_nmod_mpolyl_content(
    fq_nmod_mpoly_t g,
    const fq_nmod_mpoly_t A,
    slong num_vars,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success;
    slong i, j, off, shift;
    ulong old_shift, new_shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * Aexps = A->exps;
    slong Alen = A->length;
    fq_nmod_mpoly_struct * v;
    slong vlen, valloc;

    FLINT_ASSERT(g != A);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(0 < num_vars && num_vars < ctx->minfo->nvars);

    mpoly_gen_offset_shift_sp(&off, &shift, num_vars - 1, A->bits, ctx->minfo);

    i = 0;
    old_shift = (Aexps + N*i)[off] >> shift;

    valloc = 4;
    v = FLINT_ARRAY_ALLOC(valloc, fq_nmod_mpoly_struct);
    vlen = 0;

    v[vlen].bits = A->bits;
    v[vlen].coeffs = A->coeffs + d*i;
    v[vlen].exps = Aexps + N*i;
    v[vlen].length = i;
    v[vlen].coeffs_alloc = d*v[vlen].length;
    v[vlen].exps_alloc = N*v[vlen].length;
    vlen++;

    for (i = 1; i < Alen; old_shift = new_shift, i++)
    {
        new_shift = (Aexps + N*i)[off] >> shift;

        if (new_shift != old_shift)
            goto new_one;

        for (j = off + 1; j < N; j++)
            if ((Aexps + N*(i - 1))[j] != (Aexps + N*i)[j])
                goto new_one;

        continue;

new_one:

        v[vlen - 1].length = i - v[vlen - 1].length;
        FLINT_ASSERT(v[vlen - 1].length > 0);
        v[vlen - 1].coeffs_alloc = d*v[vlen - 1].length;
        v[vlen - 1].exps_alloc   = N*v[vlen - 1].length;

        if (vlen + 1 > valloc)
        {
            valloc += 2 + valloc/2;
            v = (fq_nmod_mpoly_struct *) flint_realloc(v, valloc*
                                                 sizeof(fq_nmod_mpoly_struct));
        }
        v[vlen].bits = A->bits;
        v[vlen].coeffs = A->coeffs + d*i;
        v[vlen].exps = Aexps + N*i;
        v[vlen].length = i;
        vlen++;
    }

    v[vlen - 1].length = i - v[vlen - 1].length;
    FLINT_ASSERT(v[vlen - 1].length > 0);
    v[vlen - 1].coeffs_alloc = d*v[vlen - 1].length;
    v[vlen - 1].exps_alloc   = N*v[vlen - 1].length;

    success = _fq_nmod_mpoly_vec_content_mpoly(g, v, vlen, ctx);

    if (success)
    {
        /* remove gen(0) ... gen(num_vars-1) from the answer */
        ulong * gexps;
        ulong mask;

        fq_nmod_mpoly_repack_bits_inplace(g, A->bits, ctx);
        gexps = g->exps;

        mask = (shift > 0) ? ((-UWORD(1)) >> (FLINT_BITS - shift)) : 0;
        for (i = 0; i < g->length; i++)
        {
            (gexps + N*i)[off] &= mask;
            for (j = off + 1; j < N; j++)
                (gexps + N*i)[j] = 0;
        }
    }

    flint_free(v);

    return success;
}

int fq_nmod_mpoly_content_vars(
    fq_nmod_mpoly_t g,
    const fq_nmod_mpoly_t A,
    slong * vars, slong num_vars,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    fq_nmod_mpolyv_t v, w;
    fq_nmod_mpoly_univar_t u;

    if (num_vars < 1)
    {
        fq_nmod_mpoly_set(g, A, ctx);
        return 1;
    }

    for (i = 0; i < num_vars; i++)
    {
        if (vars[i] >= (ulong) ctx->minfo->nvars)
            flint_throw(FLINT_ERROR, "fq_nmod_mpoly_content_vars: variable out of range");
    }

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        fq_nmod_mpoly_zero(g, ctx);
        return 1;
    }

    if (A->bits <= FLINT_BITS &&
        ctx->minfo->ord == ORD_LEX &&
        num_vars < ctx->minfo->nvars)
    {
        for (i = 0; i < num_vars; i++)
            if (vars[i] != i)
                goto do_general;

        if (g == A)
        {
            fq_nmod_mpoly_t t;
            fq_nmod_mpoly_init(t, ctx);
            success = fq_nmod_mpolyl_content(t, A, num_vars, ctx);
            fq_nmod_mpoly_swap(g, t, ctx);
            fq_nmod_mpoly_clear(t, ctx);
            return success;
        }

        return fq_nmod_mpolyl_content(g, A, num_vars, ctx);
    }

do_general:

    fq_nmod_mpolyv_init(v, ctx);

    fq_nmod_mpolyv_init(w, ctx);
    fq_nmod_mpoly_univar_init(u, ctx);

    i = 0;
    fq_nmod_mpoly_to_univar(u, A, vars[i], ctx);
    fq_nmod_mpolyv_fit_length(v, u->length, ctx);
    v->length = u->length;
    for (j = 0; j < u->length; j++)
        fq_nmod_mpoly_swap(v->coeffs + j, u->coeffs + j, ctx);

    for (i = 1; i < num_vars; i++)
    {
        w->length = 0;
        for (k = 0; k < v->length; k++)
        {
            fq_nmod_mpoly_to_univar(u, v->coeffs + k, vars[i], ctx);
            fq_nmod_mpolyv_fit_length(w, w->length + u->length, ctx);
            for (j = 0; j < u->length; j++)
            {
                fq_nmod_mpoly_swap(w->coeffs + w->length, u->coeffs + j, ctx);
                w->length++;
            }
        }
        fq_nmod_mpolyv_swap(v, w, ctx);
    }

    fq_nmod_mpoly_univar_clear(u, ctx);
    fq_nmod_mpolyv_clear(w, ctx);

    success = _fq_nmod_mpoly_vec_content_mpoly(g, v->coeffs, v->length, ctx);

    fq_nmod_mpolyv_clear(v, ctx);

    return success;
}

