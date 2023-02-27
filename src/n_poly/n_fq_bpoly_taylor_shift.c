/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


static void _n_fq_poly_taylor_shift_horner_n_fq(
    mp_limb_t * poly,
    const mp_limb_t * c,
    slong n,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, j;
    mp_limb_t * p = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            n_fq_mul(p, poly + d*(j + 1), c, ctx);
            n_fq_add(poly + d*j, poly + d*j, p, ctx);
        }
    }

    flint_free(p);
}


void n_fq_bpoly_taylor_shift_gen1_fq_nmod(
    n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_t c_,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    mp_limb_t * c = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_fq_set_fq_nmod(c, c_, ctx);
    n_fq_bpoly_set(A, B, ctx);
    for (i = A->length - 1; i >= 0; i--)
        _n_fq_poly_taylor_shift_horner_n_fq(A->coeffs[i].coeffs, c, A->coeffs[i].length, ctx);

    flint_free(c);  
}

void n_fq_bpoly_taylor_shift_gen0_fq_nmod(
    n_bpoly_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong n, i, j;
    mp_limb_t * c;
    n_poly_t t;

    if (fq_nmod_is_zero(alpha, ctx))
        return;

    c = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    n_fq_set_fq_nmod(c, alpha, ctx);

    n_poly_init(t);
    n = A->length;

    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            n_fq_poly_scalar_mul_n_fq(t, A->coeffs + j + 1, c, ctx);
            n_fq_poly_add(A->coeffs + j, A->coeffs + j, t, ctx);
        }
    }

    n_poly_clear(t);

    flint_free(c);
}



void n_fq_bpoly_taylor_shift_gen0_n_fq(
    n_fq_bpoly_t A,
    const mp_limb_t * alpha,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, j, n = A->length;
    mp_limb_t * tmp, * c, * alphainv;
    TMP_INIT;

    if (_n_fq_is_zero(alpha, d))
        return;

    TMP_START;

    tmp = (mp_limb_t *) TMP_ALLOC(d*N_FQ_MUL_INV_ITCH*sizeof(mp_limb_t));
    c = TMP_ALLOC(d*sizeof(mp_limb_t));
    alphainv = TMP_ALLOC(d*sizeof(mp_limb_t));

    _n_fq_one(c, d);
    for (i = 1; i < n; i++)
    {
        _n_fq_mul(c, c, alpha, ctx, tmp);
        if (!_n_fq_is_one(c, d))
        {
            mp_limb_t * Aic = A->coeffs[i].coeffs;
            for (j = 0; j < A->coeffs[i].length; j++)
                _n_fq_mul(Aic + d*j, Aic + d*j, c, ctx, tmp);
        }
    }

    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            n_fq_poly_add(A->coeffs + j, A->coeffs + j, A->coeffs + j + 1, ctx);
        }
    }

    _n_fq_inv(alphainv, alpha, ctx, tmp);
    _n_fq_one(c, d);
    for (i = 1; i < n; i++)
    {
        _n_fq_mul(c, c, alphainv, ctx, tmp);
        if (!_n_fq_is_one(c, d))
        {
            mp_limb_t * Aic = A->coeffs[i].coeffs;
            for (j = 0; j < A->coeffs[i].length; j++)
                _n_fq_mul(Aic + d*j, Aic + d*j, c, ctx, tmp);
        }
    }

    TMP_END;

    return;
}

