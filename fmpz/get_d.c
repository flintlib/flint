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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

#if FLINT64   /* 2^53 */
#define DOUBLE_MAX 9007199254740992L
#define DOUBLE_MIN -9007199254740992L
#else
#define DOUBLE_MAX COEFF_MAX
#define DOUBLE_MIN COEFF_MIN
#endif

extern double __gmpn_get_d(mp_limb_t *, size_t, size_t, long);

double
fmpz_get_d(const fmpz_t f)
{
    fmpz c = *f;

    if (c >= DOUBLE_MIN && c <= DOUBLE_MAX)
    {
        return (double) c;
    }
    else if (!COEFF_IS_MPZ(c))
    {
        mp_limb_t d;

        if (c > 0)
        {
            d = c;
            return __gmpn_get_d(&d, 1, 1, 0);
        }
        else
        {
            d = -c;
            return __gmpn_get_d(&d, 1, -1, 0);
        }
    }
    else
        return mpz_get_d(COEFF_TO_PTR(c));
}
