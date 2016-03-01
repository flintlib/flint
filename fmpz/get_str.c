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

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

char * fmpz_get_str(char * str, int b, const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
    {
        mpz_t z;

        flint_mpz_init_set_si(z, *f);

        if (!str) {
          str = flint_malloc(mpz_sizeinbase (z, b) + 2);
        }
        str = mpz_get_str(str, b, z);
        mpz_clear(z);
    }
    else
    {
        if (!str) {
          str = flint_malloc(mpz_sizeinbase (COEFF_TO_PTR(*f), b) + 2);
        }
        str = mpz_get_str(str, b, COEFF_TO_PTR(*f));
    }

    return str;
}

