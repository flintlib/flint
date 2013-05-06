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
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

long
fmpz_flog_ui(const fmpz_t n, ulong b)
{
    fmpz_t t;
    int sign;
    long r;

    if (fmpz_is_one(n))
        return 0;

    if (b == 2)
        return fmpz_bits(n) - 1;

    if (!COEFF_IS_MPZ(*n))
        return n_flog(*n, b);

    sign = fmpz_cmp_ui(n, b);
    if (sign <= 0)
        return (sign == 0) ? 1 : 0;

    r = fmpz_dlog(n) / log(b);

    fmpz_init(t);
    fmpz_set_ui(t, b);
    fmpz_pow_ui(t, t, r);
    sign = fmpz_cmp(t, n);

    /* Adjust down */
    if (sign > 0)
    {
        while (sign > 0)
        {
            fmpz_divexact_ui(t, t, b);
            sign = fmpz_cmp(t, n);
            r--;
        }
    }
    /* Adjust up */
    else if (sign < 0)
    {
        while (1)
        {
            fmpz_mul_ui(t, t, b);
            if (fmpz_cmp(t, n) <= 0)
                r++;
            else
                break;
        }
    }

    fmpz_clear(t);
    return r;
}
