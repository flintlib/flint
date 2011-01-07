/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 William Hart

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

long NORM(mp_srcptr R, long len)
{
    while (len != 0 && R[len - 1] == 0)
        len--;
    return len;
}

long
_nmod_poly_gcd_euclidean(mp_ptr G, mp_srcptr A, long lenA, 
                                  mp_srcptr B, long lenB, nmod_t mod)
{
    int steps = 2;
    long lenR1, lenR2 = 0, lenR3 = 0, len = 0;

    mp_ptr F, T, Q, R3 = G, R2, R1 = nmod_vec_init(FLINT_MAX(lenA - lenB + 1, lenB - 1) + 2*lenB - 3);
    R2 = R1 + lenB - 1;
    Q = R2 + lenB - 2;
    F = R1;

    _nmod_poly_divrem(Q, R1, A, lenA, B, lenB, mod);
    lenR1 = NORM(R1, lenB - 1);

    if (lenR1 > 1)
    {
        _nmod_poly_divrem(Q, R2, B, lenB, R1, lenR1, mod);
        lenR2 = NORM(R2, lenR1 - 1);
    } else
    {
        if (lenR1 == 0)
        {
            mpn_copyi(G, B, lenB);
            nmod_vec_free(F);
            return lenB;
        }
        else
        {
            G[0] = R1[0];
            nmod_vec_free(F);
            return 1;
        }
    }

    while (lenR2 > 1)
    {
        _nmod_poly_divrem(Q, R3, R1, lenR1, R2, lenR2, mod);
        lenR3 = NORM(R3, lenR2 - 1);
        if (++steps == 3) steps = 0;
        T = R1; R1 = R2; R2 = R3; R3 = T;
        lenR1 = lenR2; lenR2 = lenR3;
    }

    if (lenR2 == 1)
    {
        len = 1;
        if (steps) 
            G[0] = R2[0];
    }
    else
    {
        len = lenR1;
        if (steps != 1)
            mpn_copyi(G, R1, lenR1);
    }

    nmod_vec_free(F);
    return len;
}

void
nmod_poly_gcd_euclidean(nmod_poly_t G, 
                 const nmod_poly_t A, const nmod_poly_t B)
{
    nmod_poly_t tG;
    mp_ptr g;
    long A_len, B_len, len;

    B_len = B->length;
    A_len = A->length;
    
    if (A_len == 0)
    {
        if (B_len == 0) nmod_poly_zero(G);
        else nmod_poly_make_monic(G, B);
        return;
    } 
    else if (B_len == 0)
    {
        nmod_poly_make_monic(G, A);
        return;
    }

    if (A_len == 1 || B_len == 1)
    {
        nmod_poly_set_coeff_ui(G, 0, 1);
        G->length = 1;
        return;
    }

    if (G == A || G == B)
    {
        nmod_poly_init2(tG, A->mod.n, FLINT_MIN(A_len, B_len));
        g = tG->coeffs;
    }
    else
    {
        nmod_poly_fit_length(G, FLINT_MIN(A_len, B_len));
        g = G->coeffs;
    }

    if (A_len >= B_len)
        len = _nmod_poly_gcd_euclidean(g, A->coeffs, A_len,
                            B->coeffs, B_len, A->mod);
    else
        len = _nmod_poly_gcd_euclidean(g, B->coeffs, B_len,
                            A->coeffs, A_len, A->mod);

    if (G == A || G == B)
    {
        nmod_poly_swap(tG, G);
        nmod_poly_clear(tG);
    }
    
    G->length = len;

    if (G->length == 1)
        G->coeffs[0] = 1;
    else
        nmod_poly_make_monic(G, G);
}
