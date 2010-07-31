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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void _fmpz_poly_evaluate_horner(fmpz_t res, const fmpz * f, long len, 
                                const fmpz_t a)
{
    if (len == 0L)
        fmpz_set_ui(res, 0UL);
    else if (len == 1L)
        fmpz_set(res, f);
    else
    {
        fmpz_t t;
        const fmpz * c = f + (len - 1L);
        fmpz_init(t);
        fmpz_set(res, c);
        do {
            fmpz_mul(t, res, a);
            fmpz_add(res, --c, t);
        } while (c != f);
        fmpz_clear(t);
    }
}

void fmpz_poly_evaluate_horner(fmpz_t res, const fmpz_poly_t f, const fmpz_t a)
{
    if (res == a)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_poly_evaluate_horner(t, f->coeffs, f->length, a);
        fmpz_swap(res, t);
        fmpz_clear(t);
    }
    else
        _fmpz_poly_evaluate_horner(res, f->coeffs, f->length, a);
}

