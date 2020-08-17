/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_div(fmpz * Q, fmpz_t q, 
                    const fmpz * A, const fmpz_t a, slong lenA, 
          const fmpz * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv)
{
    slong lenQ = lenA - lenB + 1;
    ulong d;
    const fmpz * lead = B + (lenB - 1);
    
    if (lenB == 1)
    {
        _fmpq_poly_scalar_div_fmpq(Q, q, A, a, lenA, B, b);
        return;
    }
    
    /* 
       From pseudo division over Z we have 
           lead^d * A = Q * B + R
       and thus
           {A, a} = {b * Q, a * lead^d} * {B, b} + {R, a * lead^d}.
     */
    _fmpz_poly_pseudo_div(Q, &d, A, lenA, B, lenB, inv);
    
    /* 1.  lead^d == +-1.  {Q, q} = {b Q, a} up to sign */
    if (d == UWORD(0) || *lead == WORD(1) || *lead == WORD(-1))
    {
        fmpz_one(q);
        _fmpq_poly_scalar_mul_fmpz(Q, q, Q, q, lenQ, b);
        _fmpq_poly_scalar_div_fmpz(Q, q, Q, q, lenQ, a);
        
        if (*lead == WORD(-1) && d % UWORD(2))
            _fmpz_vec_neg(Q, Q, lenQ);
    }
    /* 2.  lead^d != +-1.  {Q, q} = {b Q, a lead^d} */
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
        
        fmpz_one(q);
        _fmpq_poly_scalar_mul_fmpz(Q, q, Q, q, lenQ, b);
        _fmpq_poly_scalar_div_fmpz(Q, q, Q, q, lenQ, den);
        
        fmpz_clear(den);
    }
}

void fmpq_poly_div(fmpq_poly_t Q, 
                      const fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    slong lenA, lenB, lenQ;

    if (fmpq_poly_is_zero(poly2))
    {
        flint_printf("Exception (fmpq_poly_div). Division by zero.\n");
        flint_abort();
    }

    if (poly1->length < poly2->length)
    {
        fmpq_poly_zero(Q);
        return;
    }
 
    /* Deal with aliasing */
    if (Q == poly1 || Q == poly2)
    {
        fmpq_poly_t tempQ;
        fmpq_poly_init(tempQ);
        fmpq_poly_div(tempQ, poly1, poly2);
        fmpq_poly_swap(Q, tempQ);
        fmpq_poly_clear(tempQ);
        return;
    }
    
    
    lenA = poly1->length;
    lenB = poly2->length;
    lenQ = lenA - lenB + 1;
    
    fmpq_poly_fit_length(Q, lenQ);
    
    _fmpq_poly_div(Q->coeffs, Q->den, 
                   poly1->coeffs, poly1->den, poly1->length, 
                   poly2->coeffs, poly2->den, poly2->length, NULL);
    
    _fmpq_poly_set_length(Q, lenQ);
}

