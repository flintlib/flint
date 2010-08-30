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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_signature(ulong * r1, ulong * r2, fmpz * poly, long len)
{
    fmpz *A, *B;
    long lenA, lenB;
    int s, t;
    fmpz_t g;
    
    if (len <= 2)
    {
        *r1 = (len == 2);
        *r2 = 0;
        return;
    }
    
    lenA = len;
    A = _fmpz_vec_init(lenA);
    _fmpz_poly_primitive_part(A, poly, lenA);

    lenB = lenA - 1;
	B = _fmpz_vec_init(lenB);
    _fmpz_poly_derivative(B, A, lenA);
    _fmpz_poly_primitive_part(B, B, lenB);
    
    fmpz_init(g);
    fmpz_set_ui(g, 1);
    
    s = 1;
    t = (lenA & 1L) ? -s : s;
    *r1 = 1;
    
    while (1)
	{
        int sgnA;

        _fmpz_poly_pseudo_rem_cohen(A, A, lenA, B, lenB);

        for (lenA = lenB - 1; (lenA >= 0) && (A[lenA] == 0L); lenA--) ;
        lenA++;
        
		if (lenA == 0)
		{
			printf("Exception: non-squarefree polynomial detected in fmpz_poly_signature\n");
            fmpz_clear(g);
            _fmpz_vec_clear(A, len);
            _fmpz_vec_clear(B, len - 1);
			abort();
		}
      
        if ((fmpz_sgn(B + (lenB - 1)) > 0) || ((lenA - lenB) & 1L))
            _fmpz_vec_neg(A, A, lenA);

        sgnA = fmpz_sgn(A + (lenA - 1));
		if (sgnA != s)
		{
			s = -s;
			(*r1)--;
        }
		if (sgnA != ((lenA & 1L) ? t : -t))
		{
			t = -t;
			(*r1)++;
		}

        if (lenA == 1)
        {
            *r2 = ((len - 1) - *r1) / 2;
            
            fmpz_clear(g);
            _fmpz_vec_clear(A, len);
            _fmpz_vec_clear(B, len - 1);
            return;
        }
		else
		{
            {
                fmpz * temp = A;
                A = B;
                B = temp;
            }
            {
                long temp = lenA;
                lenA = lenB;
                lenB = temp;
            }
            
            _fmpz_vec_scalar_divexact_fmpz(B, B, lenB, g);
            if (lenB > 2)
                fmpz_mul(g, A + (lenA - 1), A + (lenA - 1));
		}
	}
}

void fmpz_poly_signature(ulong * r1, ulong * r2, fmpz_poly_t poly)
{
    _fmpz_poly_signature(r1, r2, poly->coeffs, poly->length);
}

