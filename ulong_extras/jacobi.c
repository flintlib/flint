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

    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int
n_jacobi_unsigned(mp_limb_t x, mp_limb_t y)
{
    mp_limb_t a, b, temp;
    int s, exp;

    a = x;
    b = y;
    s = 1;

    if ((a < b) && (b != 1UL))
    {
        if (a == 0UL)
            return 0;

        temp = a;
        a = b;
        b = temp;

        count_trailing_zeros(exp, b);
        b >>= exp;

        /* We are only interested in values mod 8, so overflows don't matter here */
        if (((exp * (a * a - 1)) / 8) % 2 == 1UL)
            s = -s;

        /* We are only interested in values mod 4, so overflows don't matter here */
        if ((((a - 1) * (b - 1)) / 4) % 2 == 1UL)
            s = -s;
    }

    while (b != 1UL)
    {
        if ((a >> 2) < b)
        {
            temp = a - b;
            a = b;
            if (temp < b)
                b = temp;
            else if (temp < (b << 1))
                b = temp - a;
            else
                b = temp - (a << 1);
        }
        else
        {
            temp = a % b;
            a = b;
            b = temp;
        }

        if (b == 0UL)
            return 0;

        count_trailing_zeros(exp, b);
        b >>= exp;

        /* We are only interested in values mod 8, so overflows don't matter here */
        if (((exp * (a * a - 1)) / 8) % 2 == 1UL)
            s = -s;

        /* We are only interested in values mod 4, so overflows don't matter here */
        if ((((a - 1) * (b - 1)) / 4) % 2 == 1UL)
            s = -s;
    }

    return s;
}

int
n_jacobi(mp_limb_signed_t x, mp_limb_t y)
{
    if (x < 0L)
    {
        if (((y - 1) / 2) % 2 == 1UL)
            return -n_jacobi_unsigned(-x, y);
        else
            return n_jacobi_unsigned(-x, y);
    }
    else
        return n_jacobi_unsigned(x, y);
}
