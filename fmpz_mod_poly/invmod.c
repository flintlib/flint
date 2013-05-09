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
#include "fmpz_mod_poly.h"

int _fmpz_mod_poly_invmod(fmpz *A, 
                          const fmpz *B, len_t lenB, 
                          const fmpz *P, len_t lenP, const fmpz_t p)
{
    fmpz *G;
    len_t lenG;

    FMPZ_VEC_NORM(B, lenB);

    G = _fmpz_vec_init(lenB);

    lenG = _fmpz_mod_poly_gcdinv(G, A, B, lenB, P, lenP, p);

    if (lenG == 1 && !fmpz_is_one(G + 0))
    {
        fmpz_t invG;

        fmpz_init(invG);
        fmpz_invmod(invG, G + 0, p);
        _fmpz_mod_poly_scalar_mul_fmpz(A, A, lenP - 1, invG, p);
        fmpz_clear(invG);
    }

    _fmpz_vec_clear(G, lenB);

    return (lenG == 1);
}

int fmpz_mod_poly_invmod(fmpz_mod_poly_t A, 
                         const fmpz_mod_poly_t B, const fmpz_mod_poly_t P)
{
    const len_t lenB = B->length, lenP = P->length;
    fmpz *t;
    int ans;

    if (lenP < 2)
    {
        printf("Exception (fmpz_mod_poly_invmod). lenP < 2.\n");
        abort();
    }
    if (lenB == 0)
    {
        fmpz_mod_poly_zero(A);
        return 0;
    }
    if (lenB >= lenP)
    {
        fmpz_mod_poly_t T;

        fmpz_mod_poly_init(T, &B->p);
        fmpz_mod_poly_rem(T, B, P);
        ans = fmpz_mod_poly_invmod(A, T, P);
        fmpz_mod_poly_clear(T);
        return ans;
    }

    if (A != B && A != P)
    {
        fmpz_mod_poly_fit_length(A, lenP - 1);
        t = A->coeffs;
    }
    else
    {
        t = _fmpz_vec_init(lenP);
    }

    ans = _fmpz_mod_poly_invmod(t, B->coeffs, lenB, P->coeffs, lenP, &B->p);

    if (A == B || A == P)
    {
        _fmpz_vec_clear(A->coeffs, A->alloc);
        A->coeffs = t;
        A->alloc  = lenP - 1;
    }
    _fmpz_mod_poly_set_length(A, lenP - 1);
    _fmpz_mod_poly_normalise(A);
    return ans;
}

