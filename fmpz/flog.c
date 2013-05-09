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
#include "fmpz.h"
#include "ulong_extras.h"

len_t
fmpz_flog(const fmpz_t n, const fmpz_t b)
{
    fmpz_t t;
    int sign;
    len_t r;

    if (fmpz_is_one(n))
        return 0;

    if (!COEFF_IS_MPZ(*b))
        return fmpz_flog_ui(n, *b);

    sign = fmpz_cmp(n, b);
    if (sign <= 0)
        return (sign == 0) ? 1 : 0;

    r = fmpz_dlog(n) / fmpz_dlog(b);

    fmpz_init(t);
    fmpz_pow_ui(t, b, r);
    sign = fmpz_cmp(t, n);

    /* Adjust down */
    if (sign > 0)
    {
        while (sign > 0)
        {
            fmpz_divexact(t, t, b);
            sign = fmpz_cmp(t, n);
            r--;
        }
    }
    /* Adjust up */
    else if (sign < 0)
    {
        while (1)
        {
            fmpz_mul(t, t, b);
            if (fmpz_cmp(t, n) <= 0)
                r++;
            else
                break;
        }
    }

    fmpz_clear(t);
    return r;
}

