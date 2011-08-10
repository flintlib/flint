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
#include <limits.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_clog(mp_limb_t x, mp_limb_t b)
{
    if (x == 1)
    {
        return 0;
    }
    else if (x < b)
    {
        return 1;
    }
    else
    {
        mp_limb_t n;

        for (n = 1; x > b; n++)
        {
            x = (x + (b - 1)) / b;
        }

        return n;
    }
}
