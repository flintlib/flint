/*============================================================================

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

===============================================================================*/
/****************************************************************************

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart

*****************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void fmpq_poly_canonicalize(fmpq_poly_t poly)
{
    ulong length;
    
    length = poly->length;
    while (length && !poly->coeffs[length - 1]) length--;
    
    if (_fmpq_poly_den_is_one(poly))
        return;
    
    if (fmpq_poly_is_zero(poly))
    {
        fmpz_set_si(poly->den, 1);
    }
    else if (*poly->den == -1L)
    {
        _fmpz_vec_neg(poly->coeffs, poly->coeffs, poly->length);
        fmpz_set_si(poly->den, 1);
    }
}

