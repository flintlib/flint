/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void _nmod_mpoly_pow_rmul(
    nmod_mpoly_t A,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    ulong k,
    slong N,
    const ulong * cmpmask,
    nmod_t mod,
    nmod_mpoly_t T)
{
    flint_bitcnt_t bits = A->bits;

    FLINT_ASSERT(bits == T->bits);
    FLINT_ASSERT(Blen > 0);

    _nmod_mpoly_fit_length(&A->coeffs, &A->coeffs_alloc,
                           &A->exps, &A->exps_alloc, N, Blen + 2);

    if (k >= 2)
    {
        _nmod_mpoly_mul_johnson(A, Bcoeffs, Bexps, Blen,
                                   Bcoeffs, Bexps, Blen,
                                   bits, N, cmpmask, mod);
        k -= 2;
        while (k >= 1 && A->length > 0)
        {
            _nmod_mpoly_mul_johnson(T, A->coeffs, A->exps, A->length,
                                       Bcoeffs, Bexps, Blen,
                                       bits, N, cmpmask, mod);
            nmod_mpoly_swap(A, T, NULL);
            k -= 1;
        }
    }
    else if (k == 1)
    {
        FLINT_ASSERT(A->coeffs_alloc >= Blen);
        _nmod_vec_set(A->coeffs, Bcoeffs, Blen);
        mpoly_copy_monomials(A->exps, Bexps, Blen, N);
        A->length = Blen;
    }
    else
    {
        mpoly_monomial_zero(A->exps, N);
        A->coeffs[0] = 1;
        A->length = 1;
    }
}

void nmod_mpoly_pow_rmul(nmod_mpoly_t A, const nmod_mpoly_t B,
                                         ulong k, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_t T;
    nmod_mpoly_init(T, ctx);

    if (A == B)
    {
        nmod_mpoly_pow_rmul(T, A, k, ctx);
        nmod_mpoly_swap(T, A, ctx);
    }
    else
    {
        nmod_mpoly_one(A, ctx);
        while (k > 0)
        { 
            nmod_mpoly_mul_johnson(T, A, B, ctx);
            nmod_mpoly_swap(A, T, ctx);
            k -= 1;
        }
    }

    nmod_mpoly_clear(T, ctx);
}
