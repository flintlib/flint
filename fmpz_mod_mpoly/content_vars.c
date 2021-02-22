/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "fmpz_mod_mpoly_factor.h"

/*
    content wrt gen(0), ..., gen(num_vars-1)
    successful answer will be returned with g->bits == A->bits
*/
int fmpz_mod_mpolyl_content(
    fmpz_mod_mpoly_t g,
    const fmpz_mod_mpoly_t A,
    slong num_vars,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, off, shift;
    ulong old_shift, new_shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * Aexps = A->exps;
    slong Alen = A->length;
    fmpz_mod_mpoly_struct * v;
    slong vlen, valloc;

    FLINT_ASSERT(g != A);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(0 < num_vars && num_vars < ctx->minfo->nvars);

    mpoly_gen_offset_shift_sp(&off, &shift, num_vars - 1, A->bits, ctx->minfo);

    i = 0;
    old_shift = (Aexps + N*i)[off] >> shift;

    valloc = 4;
    v = FLINT_ARRAY_ALLOC(valloc, fmpz_mod_mpoly_struct);
    vlen = 0;

    v[vlen].bits = A->bits;
    v[vlen].coeffs = A->coeffs + i;
    v[vlen].exps = Aexps + N*i;
    v[vlen].coeffs_alloc = 0;
    v[vlen].exps_alloc = 0;
    v[vlen].length = i;
    v[vlen].coeffs_alloc = v[vlen].length;
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
        v[vlen - 1].coeffs_alloc = v[vlen - 1].length;
        v[vlen - 1].exps_alloc   = N*v[vlen - 1].length;

        if (vlen + 1 > valloc)
        {
            valloc += 2 + valloc/2;
            v = FLINT_ARRAY_REALLOC(v, valloc, fmpz_mod_mpoly_struct);
        }
        v[vlen].bits = A->bits;
        v[vlen].coeffs = A->coeffs + i;
        v[vlen].exps = Aexps + N*i;
        v[vlen].coeffs_alloc = 0;
        v[vlen].exps_alloc = 0;
        v[vlen].length = i;
        vlen++;
    }

    v[vlen - 1].length = i - v[vlen - 1].length;
    FLINT_ASSERT(v[vlen - 1].length > 0);
    v[vlen - 1].coeffs_alloc = v[vlen - 1].length;
    v[vlen - 1].exps_alloc   = N*v[vlen - 1].length;

    success = _fmpz_mod_mpoly_vec_content_mpoly(g, v, vlen, ctx);

    if (success)
    {
        /* remove gen(0) ... gen(num_vars-1) from the answer */
        ulong * gexps;
        ulong mask;

        fmpz_mod_mpoly_repack_bits_inplace(g, A->bits, ctx);
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

int fmpz_mod_mpoly_content_vars(
    fmpz_mod_mpoly_t g,
    const fmpz_mod_mpoly_t A,
    slong * vars, slong num_vars,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    fmpz_mod_mpolyv_t v, w;
    fmpz_mod_mpoly_univar_t u;

    if (num_vars < 1)
    {
        fmpz_mod_mpoly_set(g, A, ctx);
        return 1;
    }

    for (i = 0; i < num_vars; i++)
    {
        if (vars[i] >= (ulong) ctx->minfo->nvars)
            flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_content_vars: variable out of range");
    }

    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        fmpz_mod_mpoly_zero(g, ctx);
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
            fmpz_mod_mpoly_t t;
            fmpz_mod_mpoly_init(t, ctx);
            success = fmpz_mod_mpolyl_content(t, A, num_vars, ctx);
            fmpz_mod_mpoly_swap(g, t, ctx);
            fmpz_mod_mpoly_clear(t, ctx);
            return success;
        }

        return fmpz_mod_mpolyl_content(g, A, num_vars, ctx);
    }

do_general:

    fmpz_mod_mpolyv_init(v, ctx);

    fmpz_mod_mpolyv_init(w, ctx);
    fmpz_mod_mpoly_univar_init(u, ctx);

    i = 0;
    fmpz_mod_mpoly_to_univar(u, A, vars[i], ctx);
    fmpz_mod_mpolyv_fit_length(v, u->length, ctx);
    v->length = u->length;
    for (j = 0; j < u->length; j++)
        fmpz_mod_mpoly_swap(v->coeffs + j, u->coeffs + j, ctx);

    for (i = 1; i < num_vars; i++)
    {
        w->length = 0;
        for (k = 0; k < v->length; k++)
        {
            fmpz_mod_mpoly_to_univar(u, v->coeffs + k, vars[i], ctx);
            fmpz_mod_mpolyv_fit_length(w, w->length + u->length, ctx);
            for (j = 0; j < u->length; j++)
            {
                fmpz_mod_mpoly_swap(w->coeffs + w->length, u->coeffs + j, ctx);
                w->length++;
            }
        }
        fmpz_mod_mpolyv_swap(v, w, ctx);
    }

    fmpz_mod_mpoly_univar_clear(u, ctx);
    fmpz_mod_mpolyv_clear(w, ctx);

    success = _fmpz_mod_mpoly_vec_content_mpoly(g, v->coeffs, v->length, ctx);

    fmpz_mod_mpolyv_clear(v, ctx);

    return success;
}
