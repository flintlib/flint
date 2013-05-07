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
/*****************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
fmpz_cmp_ui(const fmpz_t f, ulong g)
{
    fmpz c = *f;

    if (!COEFF_IS_MPZ(c))    /* f is small */
    {
        if (c < 0L || g > COEFF_MAX)
            return -1;
        else 
            return c < (long) g ? -1 : c > (long) g;
    }
    else                     /* f is large */
        return mpz_cmp_ui(COEFF_TO_PTR(c), g);
}
