/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


void _n_fq_poly_mullow_(
    mp_limb_t * rop,
    const mp_limb_t * op1, slong len1,
    const mp_limb_t * op2, slong len2,
    slong n,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);
    const slong fqlen = ctx->modulus->length - 1;
    const slong pfqlen = 2*fqlen - 1;
    const nmod_t mod = ctx->mod;
    const slong rlen = len1 + len2 - 1;
    const slong m = FLINT_MIN(n, rlen);
    const slong cmlen = pfqlen*m;
    const slong clen1 = pfqlen*len1;
    const slong clen2 = pfqlen*len2;
    slong i;
    mp_limb_t * tmp;
    mp_ptr cop1, cop2, crop;

    if (len1 < 1 || len2 < 1)
    {
        _nmod_vec_zero(rop, d*n);
        return;
    }

    n_poly_stack_fit_request(St, 4);

    tmp = n_poly_stack_vec_init(St, N_FQ_REDUCE_ITCH*d);

    cop1 = n_poly_stack_vec_init(St, clen1);
    for (i = 0; i < len1; i++)
    {
        flint_mpn_copyi(cop1 + pfqlen*i, op1 + d*i, d);
        flint_mpn_zero(cop1 + pfqlen*i + d, pfqlen - d);
    }

    cop2 = n_poly_stack_vec_init(St, clen2);
    for (i = 0; i < len2; i++)
    {
        flint_mpn_copyi(cop2 + pfqlen*i, op2 + d*i, d);
        flint_mpn_zero(cop2 + pfqlen*i + d, pfqlen - d);
    }

    crop = n_poly_stack_vec_init(St, cmlen);
    if (clen1 >= clen2)
        _nmod_poly_mullow(crop, cop1, clen1, cop2, clen2, cmlen, mod);
    else
        _nmod_poly_mullow(crop, cop2, clen2, cop1, clen1, cmlen, mod);

    for (i = 0; i < m; i++)
        _n_fq_reduce2(rop + d*i, crop + pfqlen*i, ctx, tmp);

    for ( ; i < n; i++)
        _n_fq_zero(rop + d*i, d);

    n_poly_stack_vec_clear(St);
    n_poly_stack_vec_clear(St);
    n_poly_stack_vec_clear(St);
    n_poly_stack_vec_clear(St);
}


void n_fq_poly_mullow_(
    n_fq_poly_t A,
    const n_fq_poly_t B,
    const n_fq_poly_t C,
    slong order,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Blen = B->length;
    slong Clen = C->length;
    const slong m = FLINT_MIN(order, Blen + Clen - 1);
#if FLINT_WANT_ASSERT
    fq_nmod_poly_t AA, BB, CC;
#endif

    FLINT_ASSERT(n_fq_poly_is_canonical(B, ctx));
    FLINT_ASSERT(n_fq_poly_is_canonical(C, ctx));

    if (Blen < 1 || Clen < 1 || order < 1)
    {
        A->length = 0;
        return;
    }

    if (A == B || A == C)
    {
        n_fq_poly_t T;
        n_fq_poly_init(T);
        n_fq_poly_mullow_(T, B, C, order, ctx, St);
        n_fq_poly_swap(A, T);
        n_fq_poly_clear(T);
        return;
    }

#if FLINT_WANT_ASSERT
    fq_nmod_poly_init(AA, ctx);
    fq_nmod_poly_init(BB, ctx);
    fq_nmod_poly_init(CC, ctx);
    n_fq_poly_get_fq_nmod_poly(BB, B, ctx);
    n_fq_poly_get_fq_nmod_poly(CC, C, ctx);
    fq_nmod_poly_mullow(AA, BB, CC, order, ctx);
#endif

    n_poly_fit_length(A, d*m);
    _n_fq_poly_mullow_(A->coeffs, B->coeffs, Blen, C->coeffs, Clen, m, ctx, St);
    A->length = m;
    _n_fq_poly_normalise(A, d);

#if FLINT_WANT_ASSERT
    n_fq_poly_get_fq_nmod_poly(BB, A, ctx);
    FLINT_ASSERT(fq_nmod_poly_equal(BB, AA, ctx));
    fq_nmod_poly_clear(AA, ctx);
    fq_nmod_poly_clear(BB, ctx);
    fq_nmod_poly_clear(CC, ctx);
#endif

    FLINT_ASSERT(n_fq_poly_is_canonical(A, ctx));
}


void n_fq_poly_mullow(
    n_fq_poly_t A,
    const n_fq_poly_t B,
    const n_fq_poly_t C,
    slong order,
    const fq_nmod_ctx_t ctx)
{
    n_poly_stack_t St;
    n_poly_stack_init(St);
    FLINT_ASSERT(n_fq_poly_is_canonical(B, ctx));
    FLINT_ASSERT(n_fq_poly_is_canonical(C, ctx));
    n_fq_poly_mullow_(A, B, C, order, ctx, St);
    n_poly_stack_clear(St);
}
