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

    Copyright (C) 2010, 2011 Sebastian Pancratz
   
******************************************************************************/

#include "fmpq_poly.h"

#include "fmpz_poly_q.h"

void fmpz_poly_q_sub_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    if (rop == op)
    {
        fmpz_poly_q_zero(rop);
        return;
    }
    
    fmpz_poly_q_neg(rop, rop);
    fmpz_poly_q_add_in_place(rop, op);
    fmpz_poly_q_neg(rop, rop);
}

void 
fmpz_poly_q_sub(fmpz_poly_q_t rop, 
                const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
{
    fmpz_poly_t d, r2, s2;
    
    if (fmpz_poly_is_zero(op1->num))
    {
        fmpz_poly_q_neg(rop, op2);
        return;
    }
    if (fmpz_poly_is_zero(op2->num))
    {
        fmpz_poly_q_set(rop, op1);
        return;
    }
    
    if (op1 == op2)
    {
        fmpz_poly_q_zero(rop);
        return;
    }
    if (rop == op1)
    {
        fmpz_poly_q_sub_in_place(rop, op2);
        return;
    }
    if (rop == op2)
    {
        fmpz_poly_q_sub_in_place(rop, op1);
        fmpz_poly_q_neg(rop, rop);
        return;
    }
   
    /*
        From here on, we know that rop, op1 and op2 refer to distinct objects 
        in memory, although as rational functions they may still be equal
    
        XXX:  Do not maintain the remaining part of the function separately!!!
              Instead, note that this is very similar to the corresponding part 
              of the summation code.
     */

    /* Polynomials? */
    if (fmpz_poly_length(op1->den) == 1 && fmpz_poly_length(op2->den) == 1)
    {
        const len_t len1 = fmpz_poly_length(op1->num);
        const len_t len2 = fmpz_poly_length(op2->num);

        fmpz_poly_fit_length(rop->num, FLINT_MAX(len1, len2));
        _fmpq_poly_sub(rop->num->coeffs, rop->den->coeffs, 
                       op1->num->coeffs, op1->den->coeffs, len1, 
                       op2->num->coeffs, op2->den->coeffs, len2);
        _fmpz_poly_set_length(rop->num, FLINT_MAX(len1, len2));
        _fmpz_poly_set_length(rop->den, 1);
        _fmpz_poly_normalise(rop->num);
        return;
    }

    /* Denominators equal to one? */
    if (fmpz_poly_is_one(op1->den))
    {
        fmpz_poly_mul(rop->num, op1->num, op2->den);
        fmpz_poly_sub(rop->num, rop->num, op2->num);
        fmpz_poly_set(rop->den, op2->den);
        return;
    }
    if (fmpz_poly_is_one(op2->den))
    {
        fmpz_poly_mul(rop->num, op2->num, op1->den);
        fmpz_poly_sub(rop->num, op1->num, rop->num);
        fmpz_poly_set(rop->den, op1->den);
        return;
    }
    
    /* Henrici's algorithm for summation in quotient fields */

    /*
        We begin by using rop->num as a temporary variable for the gcd of the 
        two denominators' greatest common divisor
     */
    fmpz_poly_gcd(rop->num, op1->den, op2->den);
    
    if (fmpz_poly_is_one(rop->num))
    {
        fmpz_poly_mul(rop->num, op1->num, op2->den);
        fmpz_poly_mul(rop->den, op1->den, op2->num);  /* Using rop->den as temp */
        fmpz_poly_sub(rop->num, rop->num, rop->den);
        fmpz_poly_mul(rop->den, op1->den, op2->den);
    }
    else
    {
        /*
            We now copy rop->num into a new variable d, so we no longer need 
            rop->num as a temporary variable
         */
        fmpz_poly_init(d);
        fmpz_poly_swap(d, rop->num);
        
        fmpz_poly_init(r2);
        fmpz_poly_init(s2);
        
        fmpz_poly_div(r2, op1->den, d);  /* +ve leading coeff */
        fmpz_poly_div(s2, op2->den, d);  /* +ve leading coeff */
        
        fmpz_poly_mul(rop->num, op1->num, s2);
        fmpz_poly_mul(rop->den, op2->num, r2);  /* Using rop->den as temp */
        fmpz_poly_sub(rop->num, rop->num, rop->den);
        
        if (fmpz_poly_degree(rop->num) < 0)
        {
            fmpz_poly_zero(rop->den);
            fmpz_poly_set_coeff_si(rop->den, 0, 1);
        }
        else
        {
            fmpz_poly_mul(rop->den, op1->den, s2);
            
            fmpz_poly_gcd(r2, rop->num, d);
            
            if (!fmpz_poly_is_one(r2))
            {
                fmpz_poly_div(rop->num, rop->num, r2);
                fmpz_poly_div(rop->den, rop->den, r2);
            }
        }
        
        fmpz_poly_clear(d);
        fmpz_poly_clear(r2);
        fmpz_poly_clear(s2);
    }
}
