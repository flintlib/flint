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
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <limits.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

mp_limb_t n_randtest(void)
{
    ulong m;

    m = n_randlimb();

    if (m & 7UL)
    {
        return n_randbits(n_randint(FLINT_BITS + 1));
    }
    else
    {
        m >>= 3;

        if      (m % 4 == 0) return 0;
        else if (m % 4 == 1) return 1;
        else if (m % 4 == 2) return COEFF_MAX;
        else                 return LONG_MAX;
    }
}

mp_limb_t n_randtest_not_zero(void)
{
    mp_limb_t n;

    while ((n = n_randtest()) == 0) ;
    return n;
}
