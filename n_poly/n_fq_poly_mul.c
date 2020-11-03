/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


void _n_fq_poly_mul_(
    mp_limb_t * A,  /* length d*(Blen + Clen - 1) */
    const mp_limb_t * B, slong Blen,
    const mp_limb_t * C, slong Clen,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Alen = Blen + Clen - 1;
    mp_limb_t * tmp, * u;
    slong i, j;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);
    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_poly_stack_fit_request(St, 4);

    tmp = n_poly_stack_vec_init(St, 16*d);
    u = tmp + 6*d;

    if (Blen < 4 || Clen < 4)
    {
        for (i = 0; i < Alen; i++)
        {
            slong j_start = FLINT_MAX(0, i - Clen + 1);
            slong j_stop = FLINT_MIN(i, Blen - 1);
            j = j_start;
                
            if (j_stop <= j_start)
            {
                FLINT_ASSERT(j_stop == j_start);
                _n_fq_mul(A + d*i, B + d*j, C + d*(i - j), ctx, tmp);
            }
            else
            {
                _n_fq_mul2(u, B + d*j, C + d*(i - j), ctx);
                for (j++; j <= j_stop; j++)
                    _n_fq_madd2(u, B + d*j, C + d*(i - j), ctx, u + 2*d);
                _n_fq_reduce2(A + d*i, u, ctx, tmp);
            }
        }
    }
    else
    {
        mp_limb_t * xA, * xB, * xC;
        slong xAlen, xBlen, xClen, start;

        xBlen = (2*d-1)*Blen;
        xB = n_poly_stack_vec_init(St, xBlen);
        for (i = 0; i < Blen; i++)
        {
            _nmod_vec_set(xB + (2*d-1)*i, B + d*i, d);
            _nmod_vec_zero(xB + (2*d-1)*i + d, d-1);
        }
        while (xBlen > 0 && xB[xBlen - 1] == 0)
            xBlen -= 1;

        xClen = (2*d-1)*Clen;
        xC = n_poly_stack_vec_init(St, xClen);
        for (i = 0; i < Clen; i++)
        {
            _nmod_vec_set(xC + (2*d-1)*i, C + d*i, d);
            _nmod_vec_zero(xC + (2*d-1)*i + d, d-1);
        }
        while (xClen > 0 && xC[xClen - 1] == 0)
            xClen -= 1;

        if (xClen < 1 || xBlen < 1)
        {
            xAlen = 0;
            xA = n_poly_stack_vec_init(St, d); /* dummy */
        }
        else
        {
            xAlen = xClen + xBlen - 1;
            xA = n_poly_stack_vec_init(St, xAlen);
            if (xBlen >= xClen)
                _nmod_poly_mul(xA, xB, xBlen, xC, xClen, ctx->mod);
            else
                _nmod_poly_mul(xA, xC, xClen, xB, xBlen, ctx->mod);
        }

        start = 0;
        for (i = 0; i < Alen; i++, start += 2*d - 1)
        {
            if (xAlen >= start + (2*d - 1))
                _n_fq_reduce2(A + d*i, xA + start, ctx, tmp);
            else if (start < xAlen)
                _n_fq_reduce(A + d*i, xA + start, xAlen - start, ctx, tmp);
            else
                _n_fq_zero(A + d*i, d);
        }

        n_poly_stack_vec_clear(St);
        n_poly_stack_vec_clear(St);
        n_poly_stack_vec_clear(St);
    }

    n_poly_stack_vec_clear(St);
}


void n_fq_poly_mul_(
    n_fq_poly_t A,
    const n_fq_poly_t B,
    const n_fq_poly_t C,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Blen = B->length;
    slong Clen = C->length;
#if FLINT_WANT_ASSERT
    fq_nmod_poly_t AA, BB, CC;
#endif

    FLINT_ASSERT(n_fq_poly_is_canonical(B, ctx));
    FLINT_ASSERT(n_fq_poly_is_canonical(C, ctx));

    if (Blen < 1 || Clen < 1)
    {
        A->length = 0;
        return;
    }

    if (A == B || A == C)
    {
        n_fq_poly_t T;
        n_fq_poly_init(T);
        n_fq_poly_mul_(T, B, C, ctx, St);
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
    fq_nmod_poly_mul(AA, BB, CC, ctx);
#endif

    n_poly_fit_length(A, d*(Blen + Clen - 1));
    _n_fq_poly_mul_(A->coeffs, B->coeffs, Blen, C->coeffs, Clen, ctx, St);
    A->length = Blen + Clen - 1;
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


void n_fq_poly_mul(
    n_fq_poly_t A,
    const n_fq_poly_t B,
    const n_fq_poly_t C,
    const fq_nmod_ctx_t ctx)
{
    n_poly_stack_t St;
    n_poly_stack_init(St);
    FLINT_ASSERT(n_fq_poly_is_canonical(B, ctx));
    FLINT_ASSERT(n_fq_poly_is_canonical(C, ctx));
    n_fq_poly_mul_(A, B, C, ctx, St);
    n_poly_stack_clear(St);
}
