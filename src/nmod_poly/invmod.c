/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

int _nmod_poly_invmod(mp_limb_t *A,
                      const mp_limb_t *B, slong lenB,
                      const mp_limb_t *P, slong lenP, const nmod_t mod)
{
    mp_limb_t *G;
    slong lenG;

    NMOD_VEC_NORM(B, lenB);

    G = _nmod_vec_init(lenB);

    lenG = _nmod_poly_gcdinv(G, A, B, lenB, P, lenP, mod);

    if (lenG == 1 && G[0] != WORD(1))
    {
        mp_limb_t invG;

        invG = n_invmod(G[0], mod.n);
        _nmod_vec_scalar_mul_nmod(A, A, lenP - 1, invG, mod);
    }

    _nmod_vec_clear(G);

    return (lenG == 1);
}

int nmod_poly_invmod(nmod_poly_t A,
                     const nmod_poly_t B, const nmod_poly_t P)
{
    const slong lenB = B->length, lenP = P->length;
    mp_limb_t *t;
    int ans;

    if (lenP < 2)
        flint_throw(FLINT_ERROR, "lenP < 2 in %s\n", __func__);

    if (lenB == 0)
    {
        nmod_poly_zero(A);
        return 0;
    }
    if (lenB >= lenP)
    {
        nmod_poly_t T;

        nmod_poly_init(T, A->mod.n);
        nmod_poly_rem(T, B, P);
        ans = nmod_poly_invmod(A, T, P);
        nmod_poly_clear(T);
        return ans;
    }

    if (A != B && A != P)
    {
        nmod_poly_fit_length(A, lenP - 1);
        t = A->coeffs;
    }
    else
    {
        t = _nmod_vec_init(lenP);
    }

    ans = _nmod_poly_invmod(t, B->coeffs, lenB, P->coeffs, lenP, A->mod);

    if (A == B || A == P)
    {
        _nmod_vec_clear(A->coeffs);
        A->coeffs = t;
        A->alloc  = lenP - 1;
    }
    _nmod_poly_set_length(A, lenP - 1);
    _nmod_poly_normalise(A);
    return ans;
}

