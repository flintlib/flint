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
#include "fmpz_poly.h"

void
_fmpz_poly_pow_binomial(fmpz * res, const fmpz * poly, ulong e)
{
    ulong i, f;
    fmpz_t a, b, c;
    
    *a = 1L;
    *b = 1L;
    *c = 1L;
    
    fmpz_one(res);
    fmpz_one(res + e);
    
    for (i = 1UL, f = e - 1UL; i <= (e - 1UL) >> 1; i++, f--)
    {
        fmpz_mul(a, a, poly);
        fmpz_mul(b, b, poly + 1);
        fmpz_mul_ui(c, c, f + 1UL);
        fmpz_divexact_ui(c, c, i);
        
        fmpz_mul(res + i, b, c);
        fmpz_mul(res + f, a, c);
    }
    
    if ((e & 1UL) == 0UL)
    {
        fmpz_mul(a, a, poly);
        fmpz_mul(b, b, poly + 1);
        fmpz_mul_ui(c, c, f + 1UL);
        fmpz_divexact_ui(c, c, i);
        
        fmpz_mul(res + i, b, c);
        fmpz_mul(res + i, res + i, a);
        i++, f--;
    }
    
    for ( ; i <= e; i++, f--)
    {
        fmpz_mul(a, a, poly);
        fmpz_mul(b, b, poly + 1);
        
        fmpz_mul(res + i, res + i, b);
        fmpz_mul(res + f, res + f, a);
    }
    
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
}

void
fmpz_poly_pow_binomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)
{
    const len_t len = poly->length;
    len_t rlen;

    if (len != 2)
    {
        printf("Exception (fmpz_poly_pow_binomial). poly->length not equal to 2.\n");
        abort();
    }

    if (e < 3UL)
    {
        if (e == 0UL)
            fmpz_poly_set_ui(res, 1UL);
        else if (e == 1UL)
            fmpz_poly_set(res, poly);
        else  /* e == 2UL */
            fmpz_poly_sqr(res, poly);
        return;
    }
    
    rlen = (len_t) e + 1;

    if (res != poly)
    {
        fmpz_poly_fit_length(res, rlen);
        _fmpz_poly_set_length(res, rlen);
        _fmpz_poly_pow_binomial(res->coeffs, poly->coeffs, e);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        _fmpz_poly_set_length(t, rlen);
        _fmpz_poly_pow_binomial(t->coeffs, poly->coeffs, e);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
