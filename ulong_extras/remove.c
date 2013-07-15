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

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int
n_remove(mp_limb_t * n, mp_limb_t p)
{
    int exp, i;
    mp_limb_t powp[6];
    mp_limb_t quot, rem;

    if (p == 2)
    {
        count_trailing_zeros(exp, *n);
        if (exp)
            (*n) >>= exp;

        return exp;
    }

    powp[0] = p;

    for (i = 0; ; i++)
    {
        if ((*n) < powp[i])
            break;
        quot = (*n) / powp[i];
        rem = (*n) - quot * powp[i];
        if (rem != 0UL)
            break;
        powp[i + 1] = powp[i] * powp[i];
        (*n) = quot;
    }

    exp = (1 << i) - 1;

    while (i > 0)
    {
        i--;
        if ((*n) < powp[i])
            continue;
        quot = (*n) / powp[i];
        rem = (*n) - quot * powp[i];
        if (rem == 0UL)
        {
            exp += (1UL << i);
            (*n) = quot;
        }
    }

    return exp;
}
