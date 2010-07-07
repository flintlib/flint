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
#include "fmpq_poly.h"

void fmpq_poly_scalar_div_mpq(fmpq_poly_t rop, const fmpq_poly_t op, const mpq_t c)
{
    if (mpq_sgn(c) == 0)
    {
        printf("ERROR (fmpq_poly_scalar_div_mpq).  Division by 0.\n");
        abort();
    }
    
    if (mpz_cmpabs_ui(mpq_numref(c), 1UL) == 0)
    {
        fmpq_poly_scalar_mul_mpz(rop, op, mpq_denref(c));
        if (mpz_cmp_si(mpq_numref(c), -1) == 0)
            fmpq_poly_neg(rop, rop);
        return;
    }
    
    fmpq_poly_fit_length(rop, op->length);
    
    /* First set {num,den} = {denref(c),numref(c)} with the right sign. */
    /* Then proceed as in the multiplication code.                      */
    fmpz_t d, num, den;
    fmpz_init(d);
    fmpz_init(num);
    fmpz_init(den);
    fmpz_set_mpz(num, mpq_denref(c));
    fmpz_set_mpz(den, mpq_numref(c));
    if (fmpz_sgn(den) < 0)
    {
        fmpz_neg(num, num);
        fmpz_neg(den, den);
    }
    
    _fmpz_vec_content(d, op->coeffs, op->length);
    fmpz_gcd(d, d, den);
    if (*d != 1)
    {
        _fmpz_vec_scalar_divexact(rop->coeffs, op->coeffs, op->length, d);
        fmpz_divexact(den, den, d);
    }
    fmpz_gcd(d, op->den, num);
    if (*d != 1)
    {
        _fmpz_vec_scalar_mul_fmpz(rop->coeffs, op->coeffs, op->length, num);
        fmpz_mul(rop->den, op->den, den);
    }
    else
    {
        fmpz_divexact(num, num, d);
        _fmpz_vec_scalar_mul_fmpz(rop->coeffs, op->coeffs, op->length, num);
        fmpz_divexact(rop->den, op->den, den);
        fmpz_mul(rop->den, rop->den, den);
    }
    fmpz_clear(d);
    fmpz_clear(num);
    fmpz_clear(den);
    
    _fmpq_poly_set_length(rop, op->length);
    _fmpq_poly_normalise(rop);
}

