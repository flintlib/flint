/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "gr_poly.h"

static const short nmod_poly_hgcd_iter_recursive_cutoff_tab[64] =
{
    68, 88, 60, 76, 64, 64, 76, 68, 72, 68, 80, 88, 96, 108, 68, 140, 116, 108,
    96, 92, 104, 132, 140, 136, 88, 104, 120, 112, 136, 200, 188, 176, 172, 112,
    144, 156, 132, 172, 156, 148, 180, 200, 188, 164, 160, 168, 168, 168, 168,
    244, 180, 180, 200, 204, 192, 228, 212, 208, 168, 196, 188, 216, 208, 76,
};

static const short nmod_poly_gcd_hgcd_cutoff_tab[64] =
{
    493, 569, 542, 626, 626, 657, 657, 723, 759, 796, 796, 835, 835, 919, 919,
    964, 835, 796, 876, 876, 1491, 1420, 1420, 1420, 1420, 1420, 1565, 1643,
    1643, 1725, 1811, 1725, 470, 569, 517, 517, 569, 542, 569, 493, 569, 569,
    759, 542, 723, 796, 964, 835, 964, 1012, 919, 919, 964, 1062, 1062, 964,
    1012, 1115, 1012, 1115, 1115, 1170, 1228, 1643,
};

slong nmod_poly_hgcd_iter_recursive_cutoff(nmod_t mod)
{
    return nmod_poly_hgcd_iter_recursive_cutoff_tab[NMOD_BITS(mod) - 1];
}

slong nmod_poly_gcd_hgcd_cutoff(nmod_t mod)
{
    return nmod_poly_gcd_hgcd_cutoff_tab[NMOD_BITS(mod) - 1];
}


slong _nmod_poly_gcd(nn_ptr G, nn_srcptr A, slong lenA,
                              nn_srcptr B, slong lenB, nmod_t mod)
{
    if (lenB < nmod_poly_gcd_hgcd_cutoff(mod))
        return _nmod_poly_gcd_euclidean(G, A, lenA, B, lenB, mod);
    else
        return _nmod_poly_gcd_hgcd(G, A, lenA, B, lenB, mod);
}

void nmod_poly_gcd(nmod_poly_t G,
                             const nmod_poly_t A, const nmod_poly_t B)
{
    if (A->length < B->length)
    {
        nmod_poly_gcd(G, B, A);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        nmod_poly_t tG;
        nn_ptr g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            nmod_poly_zero(G);
        }
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            nmod_poly_make_monic(G, A);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                nmod_poly_init2(tG, A->mod.n, FLINT_MIN(lenA, lenB));
                g = tG->coeffs;
            }
            else
            {
                nmod_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _nmod_poly_gcd(g, A->coeffs, lenA,
                                               B->coeffs, lenB, A->mod);

            if (G == A || G == B)
            {
                nmod_poly_swap(tG, G);
                nmod_poly_clear(tG);
            }
            G->length = lenG;

            if (G->length == 1)
                G->coeffs[0] = 1;
            else
                nmod_poly_make_monic(G, G);
        }
    }
}

slong _nmod_poly_gcd_euclidean(nn_ptr G, nn_srcptr A, slong lenA,
                                        nn_srcptr B, slong lenB, nmod_t mod)
{
    slong steps;
    slong lenR1, lenR2 = 0, lenG = 0;

    nn_ptr F, R1, R2, R3 = G, T;

    if (lenB == 1)
    {
        G[0] = B[0];
        return 1;
    }

    F  = _nmod_vec_init(2*lenB - 3);
    R1 = F;
    R2 = R1 + lenB - 1;

    _nmod_poly_rem(R1, A, lenA, B, lenB, mod);
    lenR1 = lenB - 1;
    MPN_NORM(R1, lenR1);

    if (lenR1 > 1)
    {
        _nmod_poly_rem(R2, B, lenB, R1, lenR1, mod);
        lenR2 = lenR1 - 1;
        MPN_NORM(R2, lenR2);
    }
    else
    {
        if (lenR1 == 0)
        {
            flint_mpn_copyi(G, B, lenB);
            _nmod_vec_clear(F);
            return lenB;
        }
        else
        {
            G[0] = R1[0];
            _nmod_vec_clear(F);
            return 1;
        }
    }

    for (steps = 2; lenR2 > 1; steps++)
    {
        _nmod_poly_rem(R3, R1, lenR1, R2, lenR2, mod);
        lenR1 = lenR2--;
        MPN_NORM(R3, lenR2);
        T = R1; R1 = R2; R2 = R3; R3 = T;
    }

    if (lenR2 == 1)
    {
        lenG = 1;
        if (steps % 3)
            G[0] = R2[0];
    }
    else
    {
        lenG = lenR1;
        if (steps % 3 != 1)
            flint_mpn_copyi(G, R1, lenR1);
    }

    _nmod_vec_clear(F);
    return lenG;
}

void nmod_poly_gcd_euclidean(nmod_poly_t G,
                             const nmod_poly_t A, const nmod_poly_t B)
{
    if (A->length < B->length)
    {
        nmod_poly_gcd_euclidean(G, B, A);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        nmod_poly_t tG;
        nn_ptr g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            nmod_poly_zero(G);
        }
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            nmod_poly_make_monic(G, A);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                nmod_poly_init2(tG, A->mod.n, FLINT_MIN(lenA, lenB));
                g = tG->coeffs;
            }
            else
            {
                nmod_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _nmod_poly_gcd_euclidean(g, A->coeffs, lenA,
                                               B->coeffs, lenB, A->mod);

            if (G == A || G == B)
            {
                nmod_poly_swap(tG, G);
                nmod_poly_clear(tG);
            }
            G->length = lenG;

            if (G->length == 1)
                G->coeffs[0] = 1;
            else
                nmod_poly_make_monic(G, G);
        }
    }
}

slong _nmod_poly_gcd_hgcd(nn_ptr G, nn_srcptr A, slong lenA,
                                   nn_srcptr B, slong lenB, nmod_t mod)
{
    slong inner_cutoff = nmod_poly_hgcd_iter_recursive_cutoff(mod);
    slong outer_cutoff = nmod_poly_gcd_hgcd_cutoff(mod);
    slong lenG = 0;
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_gcd_hgcd(G, &lenG, A, lenA, B, lenB, inner_cutoff, outer_cutoff, ctx));
    return lenG;
}

void nmod_poly_gcd_hgcd(nmod_poly_t G,
                             const nmod_poly_t A, const nmod_poly_t B)
{
    if (A->length < B->length)
    {
        nmod_poly_gcd_hgcd(G, B, A);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        nmod_poly_t tG;
        nn_ptr g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            nmod_poly_zero(G);
        }
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            nmod_poly_make_monic(G, A);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                nmod_poly_init2(tG, A->mod.n, FLINT_MIN(lenA, lenB));
                g = tG->coeffs;
            }
            else
            {
                nmod_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _nmod_poly_gcd_hgcd(g, A->coeffs, lenA,
                                               B->coeffs, lenB, A->mod);

            if (G == A || G == B)
            {
                nmod_poly_swap(tG, G);
                nmod_poly_clear(tG);
            }
            G->length = lenG;

            if (G->length == 1)
                G->coeffs[0] = 1;
            else
                nmod_poly_make_monic(G, G);
        }
    }
}
