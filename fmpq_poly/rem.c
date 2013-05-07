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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_rem(fmpz * R, fmpz_t r, 
                    const fmpz * A, const fmpz_t a, long lenA, 
                    const fmpz * B, const fmpz_t b, long lenB)
{
    long lenR;
    ulong d;
    const fmpz * lead = B + (lenB - 1);
    
    if (lenB == 1)
    {
        fmpz_one(r);
        return;
    }
    
    /* 
       From pseudo division over Z we have 
           lead^d * A = Q * B + R
       and thus
           {A, a} = {b * Q, a * lead^d} * {B, b} + {R, a * lead^d}.
     */
    _fmpz_poly_pseudo_rem(R, &d, A, lenA, B, lenB);
    
    /* Determine the actual length of R */
    for (lenR = lenB - 2; lenR >= 0 && !R[lenR]; lenR--) ;
    lenR++;
    
    /* 1.  lead^d == +-1.  {R, r} = {R, a} up to sign */
    if (d == 0UL || *lead == 1L || *lead == -1L)
    {
        fmpz_one(r);
        if (lenR > 0)
            _fmpq_poly_scalar_div_fmpz(R, r, R, r, lenR, a);
        
        if (*lead == -1L && d % 2UL)
            _fmpz_vec_neg(R, R, lenR);
    }
    /* 2.  lead^d != +-1.  {R, r} = {R, a lead^d} */
    else
    {
        /*
           TODO:  Improve this.  Clearly we do not need to compute 
           den = a lead^d in many cases, but can determine the GCD from 
           lead alone already.
         */
        fmpz_t den;
        fmpz_init(den);
        fmpz_pow_ui(den, lead, d);
        fmpz_mul(den, a, den);
        
        fmpz_one(r);
        if (lenR > 0)
            _fmpq_poly_scalar_div_fmpz(R, r, R, r, lenR, den);
        
        fmpz_clear(den);
    }
}

void fmpq_poly_rem(fmpq_poly_t R, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    long lenA, lenB, lenR;

    if (fmpq_poly_is_zero(poly2))
    {
        printf("Exception (fmpq_poly_rem). Division by zero.\n");
        abort();
    }

    if (poly1->length < poly2->length)
    {
        fmpq_poly_set(R, poly1);
        return;
    }

    /* Deal with aliasing */
    if (R == poly1 || R == poly2)
    {
        fmpq_poly_t tempR;
        fmpq_poly_init(tempR);
        fmpq_poly_rem(tempR, poly1, poly2);
        fmpq_poly_swap(R, tempR);
        fmpq_poly_clear(tempR);
        return;
    }
    
    lenA = poly1->length;
    lenB = poly2->length;
    lenR = lenB - 1;
    
    fmpq_poly_fit_length(R, lenA);  /* XXX: Need at least that much space */
    
    _fmpq_poly_rem(R->coeffs, R->den, 
                   poly1->coeffs, poly1->den, poly1->length, 
                   poly2->coeffs, poly2->den, poly2->length);
    
    _fmpq_poly_set_length(R, lenR);
    _fmpq_poly_normalise(R);
}

