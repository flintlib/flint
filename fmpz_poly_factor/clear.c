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
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly_factor.h"

void fmpz_poly_factor_clear(fmpz_poly_factor_t fac)
{
    if (fac->alloc)
    {
        long i;

        for (i = 0; i < fac->alloc; i++)
        {
            fmpz_poly_clear(fac->p + i);
        }

        fmpz_clear(&(fac->c));
        flint_free(fac->p);
        flint_free(fac->exp);
        fac->p   = NULL;
        fac->exp = NULL;
    }
}

