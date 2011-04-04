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

    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_powm_ui(fmpz_t f, const fmpz_t g, ulong exp, const fmpz_t m)
{
    if (fmpz_sgn(m) <= 0)
    {
        printf("Exception (fmpz_powm_ui).  Modulus is less than 1.\n");
        abort();
    }

    if (fmpz_is_one(m))
    {
        fmpz_zero(f);
        return;
    }

    if (exp == 0)
    {
        fmpz_set_ui(f, 1);
        return;
    }

    /* TODO:  Implement this properly! */
    {
        fmpz_t copy;

        fmpz_init(copy);
        fmpz_set(copy, m);

        fmpz_pow_ui(f, g, exp);
        fmpz_mod(f, f, copy);

        fmpz_clear(copy);
    }
}

