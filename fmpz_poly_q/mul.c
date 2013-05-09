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

void fmpz_poly_q_mul(fmpz_poly_q_t rop, 
                     const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
{
    if (fmpz_poly_q_is_zero(op1) || fmpz_poly_q_is_zero(op2))
    {
        fmpz_poly_q_zero(rop);
        return;
    }
    
    if (op1 == op2)
    {
        fmpz_poly_pow(rop->num, op1->num, 2);
        fmpz_poly_pow(rop->den, op1->den, 2);
        return;
    }
    if (rop == op1 || rop == op2)
    {
        fmpz_poly_q_t t;

        fmpz_poly_q_init(t);
        fmpz_poly_q_mul(t, op1, op2);
        fmpz_poly_q_swap(rop, t);
        fmpz_poly_q_clear(t);
        return; 
    }

    /*
        From here on, we may assume that rop, op1 and op2 refer to distinct 
        objects in memory, and that op1 and op2 are non-zero
     */

    /* Polynomials? */
    if (fmpz_poly_length(op1->den) == 1 && fmpz_poly_length(op2->den) == 1)
    {
        const len_t len1 = fmpz_poly_length(op1->num);
        const len_t len2 = fmpz_poly_length(op2->num);

        fmpz_poly_fit_length(rop->num, len1 + len2 - 1);
        if (len1 >= len2)
        {
            _fmpq_poly_mul(rop->num->coeffs, rop->den->coeffs, 
                           op1->num->coeffs, op1->den->coeffs, len1, 
                           op2->num->coeffs, op2->den->coeffs, len2);
        }
        else
        {
            _fmpq_poly_mul(rop->num->coeffs, rop->den->coeffs, 
                           op2->num->coeffs, op2->den->coeffs, len2, 
                           op1->num->coeffs, op1->den->coeffs, len1);
        }
        _fmpz_poly_set_length(rop->num, len1 + len2 - 1);
        _fmpz_poly_set_length(rop->den, 1);

        return;
    }
    
    fmpz_poly_gcd(rop->num, op1->num, op2->den);
    
    if (fmpz_poly_is_one(rop->num))
    {
        fmpz_poly_gcd(rop->den, op2->num, op1->den);
        
        if (fmpz_poly_is_one(rop->den))
        {
            fmpz_poly_mul(rop->num, op1->num, op2->num);
            fmpz_poly_mul(rop->den, op1->den, op2->den);
        }
        else
        {
            fmpz_poly_div(rop->num, op2->num, rop->den);
            fmpz_poly_mul(rop->num, op1->num, rop->num);
            fmpz_poly_div(rop->den, op1->den, rop->den);
            fmpz_poly_mul(rop->den, rop->den, op2->den);
        }
    }
    else
    {
        fmpz_poly_gcd(rop->den, op2->num, op1->den);
        
        if (fmpz_poly_is_one(rop->den))
        {
            fmpz_poly_div(rop->den, op2->den, rop->num);
            fmpz_poly_mul(rop->den, op1->den, rop->den);
            fmpz_poly_div(rop->num, op1->num, rop->num);
            fmpz_poly_mul(rop->num, rop->num, op2->num);
        }
        else
        {
            fmpz_poly_t t, u;

            fmpz_poly_init(t);
            fmpz_poly_init(u);
            fmpz_poly_div(t, op1->num, rop->num);
            fmpz_poly_div(u, op2->den, rop->num);
            fmpz_poly_div(rop->num, op2->num, rop->den);
            fmpz_poly_mul(rop->num, t, rop->num);
            fmpz_poly_div(rop->den, op1->den, rop->den);
            fmpz_poly_mul(rop->den, rop->den, u);
            fmpz_poly_clear(t);
            fmpz_poly_clear(u);
        }
    }
}
