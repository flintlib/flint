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
    Copyright (C) 2010 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpq_poly.h"

void fmpq_poly_get_coeff_fmpq(fmpq_t x, const fmpq_poly_t poly, long n)
{
    if (n >= poly->length)  /* Coefficient is beyond the end of poly */
    {
        fmpq_zero(x);
        return;
    }
    
    fmpz_set(fmpq_numref(x), poly->coeffs + n);
    fmpz_set(fmpq_denref(x), poly->den);
    fmpq_canonicalise(x);
}

