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

mp_limb_t
n_invmod(mp_limb_t x, mp_limb_t y)
{
    mp_limb_signed_t v1 = UWORD(0);
    mp_limb_signed_t v2 = UWORD(1);
    mp_limb_signed_t t2;
    mp_limb_t u3, v3;
    mp_limb_t quot, rem;

    u3 = y, v3 = x;

    if (v3 > u3)
    {
        rem = u3;
        u3 = v3;
        t2 = v2;
        v2 = v1;
        v1 = t2;
        v3 = rem;
    }

    if ((mp_limb_signed_t) (y & x) < WORD(0))  /* y and x both have top bit set */
    {
        quot = u3 - v3;
        t2 = v2;
        u3 = v3;
        v2 = v1 - v2;
        v1 = t2;
        v3 = quot;
    }

    while ((mp_limb_signed_t) (v3 << 1) < WORD(0))  /* second value has second msb set */
    {
        quot = u3 - v3;
        if (quot < v3)
        {
            t2 = v2;
            u3 = v3;
            v2 = v1 - v2;
            v1 = t2;
            v3 = quot;
        }
        else if (quot < (v3 << 1))
        {
            u3 = v3;
            t2 = v2;
            v2 = v1 - (v2 << 1);
            v1 = t2;
            v3 = quot - u3;
        }
        else
        {
            u3 = v3;
            t2 = v2;
            v2 = v1 - 3 * v2;
            v1 = t2;
            v3 = quot - (u3 << 1);
        }
    }

    while (v3)
    {
        if (u3 < (v3 << 2))  /* overflow not possible due to top 2 bits of v3 not being set */
        {
            quot = u3 - v3;
            if (quot < v3)
            {
                t2 = v2;
                u3 = v3;
                v2 = v1 - v2;
                v1 = t2;
                v3 = quot;
            }
            else if (quot < (v3 << 1))
            {
                u3 = v3;
                t2 = v2;
                v2 = v1 - (v2 << 1);
                v1 = t2;
                v3 = quot - u3;
            }
            else
            {
                u3 = v3;
                t2 = v2;
                v2 = v1 - 3 * v2;
                v1 = t2;
                v3 = quot - (u3 << 1);
            }
        }
        else
        {
            quot = u3 / v3;
            rem = u3 - v3 * quot;
            u3 = v3;
            t2 = v2;
            v2 = v1 - quot * v2;
            v1 = t2;
            v3 = rem;
        }
    }

    if (v1 < WORD(0))
        v1 += y;

    return v1;
}
