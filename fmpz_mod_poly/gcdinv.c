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
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

len_t _fmpz_mod_poly_gcdinv(fmpz *G, fmpz *S, 
                           const fmpz *A, len_t lenA, const fmpz *B, len_t lenB, 
                           const fmpz_t p)
{
    fmpz *T;
    fmpz_t inv;
    len_t ans;

    T = _fmpz_vec_init(lenA - 1);
    fmpz_init(inv);
    fmpz_invmod(inv, A + (lenA - 1), p);

    ans = _fmpz_mod_poly_xgcd(G, T, S, B, lenB, A, lenA, inv, p);

    fmpz_clear(inv);
    _fmpz_vec_clear(T, lenA - 1);

    return ans;
}

void fmpz_mod_poly_gcdinv(fmpz_mod_poly_t G, fmpz_mod_poly_t S, 
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    const len_t lenA = A->length, lenB = B->length;

    if (lenB < 2)
    {
        printf("Exception (fmpz_mod_poly_gcdinv). lenB < 2.\n");
        abort();
    }
    if (lenA >= lenB)
    {
        fmpz_mod_poly_t T;

        fmpz_mod_poly_init(T, &A->p);
        fmpz_mod_poly_rem(T, A, B);
        fmpz_mod_poly_gcdinv(G, S, T, B);
        fmpz_mod_poly_clear(T);
        return;
    }

    if (lenA == 0)
    {
        fmpz_mod_poly_zero(G);
        fmpz_mod_poly_zero(S);
    } 
    else
    {
        fmpz *g, *s;
        len_t lenG;

        if (G == A || G == B)
        {
            g = _fmpz_vec_init(lenA);
        }
        else
        {
            fmpz_mod_poly_fit_length(G, lenA);
            g = G->coeffs;
        }
        if (S == A || S == B)
        {
            s = _fmpz_vec_init(lenB - 1);
        }
        else
        {
            fmpz_mod_poly_fit_length(S, lenB - 1);
            s = S->coeffs;
        }

        lenG = _fmpz_mod_poly_gcdinv(g, s, 
            A->coeffs, lenA, B->coeffs, lenB, &A->p);

        if (G == A || G == B)
        {
            _fmpz_vec_clear(G->coeffs, G->alloc);
            G->coeffs = g;
            G->alloc  = lenA;
        }
        if (S == A || S == B)
        {
            _fmpz_vec_clear(S->coeffs, S->alloc);
            S->coeffs = s;
            S->alloc  = lenB - 1;
        }

        _fmpz_mod_poly_set_length(G, lenG);
        _fmpz_mod_poly_set_length(S, lenB - lenG);
        _fmpz_mod_poly_normalise(S);

        if (!fmpz_is_one(fmpz_mod_poly_lead(G)))
        {
            fmpz_t inv;

            fmpz_init(inv);
            fmpz_invmod(inv, fmpz_mod_poly_lead(G), &A->p);
            fmpz_mod_poly_scalar_mul_fmpz(G, G, inv);
            fmpz_mod_poly_scalar_mul_fmpz(S, S, inv);
            fmpz_clear(inv);
        }
    }
}

