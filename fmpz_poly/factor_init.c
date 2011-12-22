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

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_factor_init(fmpz_poly_factor_t fac)
{
   fac->factors   = NULL;
   fac->exponents = NULL;
   fac->alloc     = 0;
   fac->length    = 0;
}

void fmpz_poly_factor_init2(fmpz_poly_factor_t fac, long alloc)
{
    if (alloc)
    {
        long i;

        fac->factors   = malloc(alloc * sizeof(fmpz_poly_struct));
        fac->exponents = malloc(alloc * sizeof(long));

        for (i = 0; i < alloc; i++)
        {
            fmpz_poly_init(fac->factors + i);
            fac->exponents[i] = 0L;
        }
    }
    else
    {
        fac->factors   = NULL;
        fac->exponents = NULL;
    }

    fac->alloc  = alloc;
    fac->length = 0;
}

