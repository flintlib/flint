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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz *poly, long len, 
                                  const fmpz_t a, const fmpz_t p)
{
    if (len == 0)
    {
        fmpz_zero(res);
    }
    else if (len == 1 || fmpz_is_zero(a))
    {
        fmpz_set(res, poly);
    }
    else
    {
        long i = len - 1;
        fmpz_t t;

        fmpz_init(t);
        fmpz_set(res, poly + i);
        for (i = len - 2; i >= 0; i--)
        {
            fmpz_mul(t, res, a);
            fmpz_mod(t, t, p);
            fmpz_add(res, poly + i, t);
        }
        fmpz_clear(t);

        if (fmpz_cmpabs(res, p) >= 0)
            fmpz_sub(res, res, p);
    }
}

void fmpz_mod_poly_evaluate_fmpz(fmpz_t res, 
                                 const fmpz_mod_poly_t poly, const fmpz_t a)
{
    if (res == a)
    {
        fmpz_t t;

        fmpz_init(t);
        _fmpz_mod_poly_evaluate_fmpz(t, poly->coeffs, poly->length, 
                                     a, &(poly->p));
        fmpz_swap(res, t);
        fmpz_clear(t);
    }
    else
    {
        _fmpz_mod_poly_evaluate_fmpz(res, poly->coeffs, poly->length, 
                                     a, &(poly->p));
    }
}

