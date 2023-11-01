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

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

slong _nmod_poly_gcdinv(mp_limb_t *G, mp_limb_t *S,
                        const mp_limb_t *A, slong lenA,
                        const mp_limb_t *B, slong lenB,
                        const nmod_t mod)
{
    mp_limb_t *T;
    slong ans;

    T = _nmod_vec_init(lenA - 1);

    ans = _nmod_poly_xgcd(G, T, S, B, lenB, A, lenA, mod);

    _nmod_vec_clear(T);

    return ans;
}

void nmod_poly_gcdinv(nmod_poly_t G, nmod_poly_t S,
                      const nmod_poly_t A, const nmod_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;

    if (lenB < 2)
        flint_throw(FLINT_ERROR, "lenB < 2 in %s\n", __func__);

    if (lenA >= lenB)
    {
        nmod_poly_t T;

        /* TODO: We can probably use init_preinv here */
        nmod_poly_init(T, A->mod.n);
        nmod_poly_rem(T, A, B);
        nmod_poly_gcdinv(G, S, T, B);
        nmod_poly_clear(T);
        return;
    }

    if (lenA == 0)
    {
        nmod_poly_zero(G);
        nmod_poly_zero(S);
    }
    else
    {
        mp_limb_t *g, *s;
        slong lenG;

        if (G == A || G == B)
        {
            g = _nmod_vec_init(lenA);
        }
        else
        {
            nmod_poly_fit_length(G, lenA);
            g = G->coeffs;
        }
        if (S == A || S == B)
        {
            s = _nmod_vec_init(lenB - 1);
        }
        else
        {
            nmod_poly_fit_length(S, lenB - 1);
            s = S->coeffs;
        }

        lenG = _nmod_poly_gcdinv(g, s,
                                 A->coeffs, lenA,
                                 B->coeffs, lenB,
                                 A->mod);

        if (G == A || G == B)
        {
            _nmod_vec_clear(G->coeffs);
            G->coeffs = g;
            G->alloc  = lenA;
        }
        if (S == A || S == B)
        {
            _nmod_vec_clear(S->coeffs);
            S->coeffs = s;
            S->alloc  = lenB - 1;
        }

        _nmod_poly_set_length(G, lenG);
        _nmod_poly_set_length(S, lenB - lenG);
        _nmod_poly_normalise(S);

        if (nmod_poly_lead(G)[0] != WORD(1))
        {
            mp_limb_t inv;

            inv = n_invmod(nmod_poly_lead(G)[0], A->mod.n);
            nmod_poly_scalar_mul_nmod(G, G, inv);
            nmod_poly_scalar_mul_nmod(S, S, inv);
        }
    }
}

