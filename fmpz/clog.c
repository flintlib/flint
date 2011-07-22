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
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

long fmpz_clog(const fmpz_t x, const fmpz_t b)
{
    if (!COEFF_IS_MPZ(*b))  /* b is small */
    {
        return fmpz_clog_ui(x, *b);
    }
    else  /* b is large */
    {
        if (fmpz_is_one(x))
        {
            return 0;
        }
        else if (fmpz_cmp(x, b) < 0)
        {
            return 1;
        }
        else
        {
            long n;
            fmpz_t t;

            fmpz_init(t);
            fmpz_set(t, b);
            for (n = 1; fmpz_cmp(t, x) < 0; n++)
            {
                fmpz_mul(t, t, b);
            }
            fmpz_clear(t);

            return n;
        }
    }
}

