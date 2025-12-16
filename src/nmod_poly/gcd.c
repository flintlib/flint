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

static const short nmod_poly_hgcd_outer_cutoff_tab[64] =
{
    475, 432, 341, 412, 697, 697, 633, 633, 664, 603, 731, 931, 767, 845, 731,
    731, 697, 697, 845, 887, 697, 1129, 931, 1025, 1185, 1025, 1185, 1244, 1371,
    1585, 1585, 1510, 375, 498, 432, 375, 453, 412, 325, 412, 498, 475, 498,
    393, 412, 310, 393, 432, 393, 575, 633, 633, 664, 548, 475, 845, 845, 548,
    767, 522, 603, 575, 731, 1510,
};

static const short nmod_poly_gcd_hgcd_cutoff_tab[64] =
{
    470, 657, 626, 689, 689, 723, 796, 876, 835, 876, 919, 964, 964, 1012, 1170,
    1228, 919, 1115, 1115, 1115, 1491, 1491, 1643, 1725, 1811, 1725, 1811, 1725,
    1901, 1901, 1901, 1170, 1115, 1170, 1228, 1170, 1228, 1170, 1228, 1289,
    1289, 1228, 1289, 1170, 1353, 1643, 1725, 1811, 1725, 1643, 1901, 1901,
    1901, 1811, 1901, 1901, 1901, 1996, 1901, 1901, 1996, 2199, 1289, 1725,
};

slong nmod_poly_hgcd_iter_recursive_cutoff(nmod_t mod)
{
    return nmod_poly_hgcd_iter_recursive_cutoff_tab[NMOD_BITS(mod) - 1];
}

slong nmod_poly_hgcd_outer_cutoff(nmod_t mod)
{
    return nmod_poly_hgcd_outer_cutoff_tab[NMOD_BITS(mod) - 1];
}

slong nmod_poly_gcd_hgcd_cutoff(nmod_t mod)
{
    return nmod_poly_gcd_hgcd_cutoff_tab[NMOD_BITS(mod) - 1];
}


slong _nmod_poly_gcd(nn_ptr G, nn_srcptr A, slong lenA,
                              nn_srcptr B, slong lenB, nmod_t mod)
{
    if (lenB < nmod_poly_gcd_hgcd_cutoff(mod))
    {
        if (NMOD_POLY_GCD_EUCLIDEAN_USE_REDC_FAST(lenB, mod))
            return _nmod_poly_gcd_euclidean_redc_fast(G, A, lenA, B, lenB, mod);
        else
            return _nmod_poly_gcd_euclidean(G, A, lenA, B, lenB, mod);
    }
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
    slong outer_cutoff = nmod_poly_hgcd_outer_cutoff(mod);
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

slong _nmod_poly_gcd_euclidean_redc_fast(nn_ptr G, nn_srcptr A, slong lenA,
                                        nn_srcptr B, slong lenB, nmod_t mod)
{
    slong lenR, lenG, i;
    nn_ptr R;
    TMP_INIT;

    FLINT_ASSERT(mod.n % 2);
    FLINT_ASSERT(mod.n != 1);
    FLINT_ASSERT(NMOD_BITS(mod) <= FLINT_BITS - 2);

    TMP_START;
    R = TMP_ALLOC((lenB - 1) * sizeof(ulong));

    /* Do initial remainder in the original representation, mainly to deal
       efficiently with unbalanced input. */
    _nmod_poly_rem(R, A, lenA, B, lenB, mod);
    lenR = lenB - 1;
    MPN_NORM(R, lenR);

    if (lenR == 0)
    {
        flint_mpn_copyi(G, B, lenB);
        lenG = lenB;
    }
    else if (lenR == 1)
    {
        G[0] = R[0];
        lenG = 1;
    }
    else
    {
        gr_ctx_t ctx;
        GR_MUST_SUCCEED(gr_ctx_init_nmod_redc_fast(ctx, mod.n));
        GR_MUST_SUCCEED(_gr_poly_gcd_euclidean(G, &lenG, B, lenB, R, lenR, ctx));
        /* The redc scaling factor does not matter for the GCD, but we need to
           convert back to normalised residues.  */
        for (i = 0; i < lenG; i++)
            G[i] = nmod_redc_fast_normalise(G[i], GR_NMOD_REDC_CTX(ctx));
    }

    TMP_END;
    return lenG;
}

void nmod_poly_gcd_euclidean_redc_fast(nmod_poly_t G,
                             const nmod_poly_t A, const nmod_poly_t B)
{
    if (A->length < B->length)
    {
        nmod_poly_gcd_euclidean_redc_fast(G, B, A);
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

            lenG = _nmod_poly_gcd_euclidean_redc_fast(g, A->coeffs, lenA,
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
