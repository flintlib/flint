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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_rescale(fmpz * res, fmpz_t denr, 
                   const fmpz * poly, const fmpz_t den, len_t len, 
                   const fmpz_t xnum, const fmpz_t xden)
{
    if (len < 2L)
    {
        if (res != poly)
        {
            _fmpz_vec_set(res, poly, len);
            fmpz_set(denr, den);
        }
    }
    else
    {
        len_t i;
        fmpz_t t;
        
        fmpz_init(t);
        
        fmpz_one(t);
        fmpz_set(res, poly);
        for (i = 1L; i < len; i++)
        {
            fmpz_mul(t, t, xnum);
            fmpz_mul(res + i, poly + i, t);
        }
        fmpz_one(t);
        for (i = len - 2L; i >= 0L; i--)
        {
            fmpz_mul(t, t, xden);
            fmpz_mul(res + i, res + i, t);
        }
        fmpz_mul(denr, den, t);
        
        fmpz_clear(t);
        
        _fmpq_poly_canonicalise(res, denr, len);
    }
}

void fmpq_poly_rescale(fmpq_poly_t res, const fmpq_poly_t poly, const fmpq_t x)
{
    if (fmpq_is_zero(x))
    {
        fmpq_poly_zero(res);
    }
    else if (poly->length < 2L)
    {
        fmpq_poly_set(res, poly);
    }
    else
    {
        fmpq_poly_fit_length(res, poly->length);
        _fmpq_poly_rescale(res->coeffs, res->den, 
                           poly->coeffs, poly->den, poly->length, 
                           fmpq_numref(x), fmpq_denref(x));
        _fmpq_poly_set_length(res, poly->length);
    }
}
