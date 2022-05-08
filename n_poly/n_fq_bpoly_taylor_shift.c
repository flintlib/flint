/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# else
#  if HAVE_ALLOCA_H
#   include <alloca.h>
#  else
#   if _MSC_VER
#    include <malloc.h>
#    define alloca _alloca
#   else
#    ifdef __DECC
#     define alloca(x) __ALLOCA(x)
#    else
#     ifdef BSD
#      include <stdlib.h>
#     else
#      error Could not find alloca
#     endif
#    endif
#   endif
#  endif
# endif
#endif

#include "flint-impl.h"
#include "n_poly.h"

static void _n_fq_poly_taylor_shift_horner_n_fq(
    ulong * poly,
    const ulong * c,
    slong n,
    const fq_nmod_ctx_t ctx)
{
    slong d = FQ_NMOD_CTX_DEGREE(ctx);
    slong i, j;
    ulong * p = FLINT_ARRAY_ALLOC(d, ulong);

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
    slong d = FQ_NMOD_CTX_DEGREE(ctx);
    slong i;
    ulong * c = FLINT_ARRAY_ALLOC(d, ulong);

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
    slong d = FQ_NMOD_CTX_DEGREE(ctx);
    slong n, i, j;
    ulong * c;
    n_poly_t t;

    if (FQ_NMOD_IS_ZERO(alpha, ctx))
        return;

    c = FLINT_ARRAY_ALLOC(d, ulong);
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
    const ulong * alpha,
    const fq_nmod_ctx_t ctx)
{
    slong d = FQ_NMOD_CTX_DEGREE(ctx);
    slong i, j, n = A->length;
    ulong * tmp, * c, * alphainv;
    TMP_INIT;

    if (_n_fq_is_zero(alpha, d))
        return;

    TMP_START;

    tmp = (ulong *) TMP_ALLOC(d*N_FQ_MUL_INV_ITCH*sizeof(ulong));
    c = TMP_ALLOC(d*sizeof(ulong));
    alphainv = TMP_ALLOC(d*sizeof(ulong));

    _n_fq_one(c, d);
    for (i = 1; i < n; i++)
    {
        _n_fq_mul(c, c, alpha, ctx, tmp);
        if (!_n_fq_is_one(c, d))
        {
            ulong * Aic = A->coeffs[i].coeffs;
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
            ulong * Aic = A->coeffs[i].coeffs;
            for (j = 0; j < A->coeffs[i].length; j++)
                _n_fq_mul(Aic + d*j, Aic + d*j, c, ctx, tmp);
        }
    }

    TMP_END;

    return;
}

