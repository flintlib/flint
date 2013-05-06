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
n_remove2_precomp(mp_limb_t * n, mp_limb_t p, double ppre)
{
    int exp = 0;
    mp_limb_t quot, rem = 0UL;

    if (p == 2)
    {
        count_trailing_zeros(exp, *n);
        if (exp)
            (*n) >>= exp;

        return exp;
    }

    do
    {
        if ((*n) < p)
            break;
        rem = n_divrem2_precomp(&quot, *n, p, ppre);
        if (rem)
            break;
        exp++;
        (*n) = quot;
    } while (1);

    return exp;
}
