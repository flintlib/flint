/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

slong _fmpz_mod_poly_gcdinv(fmpz *G, fmpz *S, 
                           const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                           const fmpz_t p)
{
    fmpz *T;
    fmpz_t inv;
    slong ans;

    fmpz_init(inv);
    fmpz_invmod(inv, A + (lenA - 1), p);

    if (lenB < 16)
	{
	    ans = _fmpz_mod_poly_gcdinv_euclidean(G, S,
		                                   A, lenA, B, lenB, inv, p);
	} else
	{
	    T = _fmpz_vec_init(lenA - 1);
    
	    ans = _fmpz_mod_poly_xgcd(G, T, S, B, lenB, A, lenA, inv, p);
		
		_fmpz_vec_clear(T, lenA - 1);
    }
	
    fmpz_clear(inv);

    return ans;
}

void fmpz_mod_poly_gcdinv(fmpz_mod_poly_t G, fmpz_mod_poly_t S, 
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;

    if (lenB < 2)
    {
        flint_printf("Exception (fmpz_mod_poly_gcdinv). lenB < 2.\n");
        flint_abort();
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
        slong lenG;

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

