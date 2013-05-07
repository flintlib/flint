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
n_powmod_ui_precomp(mp_limb_t a, mp_limb_t exp, mp_limb_t n, double npre)
{
    mp_limb_t x, y;

    if (n == 1UL)
        return 0L;

    x = 1UL;
    y = a;

    while (exp)
    {
        if (exp & 1L)
            x = n_mulmod_precomp(x, y, n, npre);
        exp >>= 1;
        if (exp)
            y = n_mulmod_precomp(y, y, n, npre);
    }

    return x;
}

mp_limb_t
n_powmod_precomp(mp_limb_t a, mp_limb_signed_t exp, mp_limb_t n, double npre)
{
    if (exp < 0)
    {
        a = n_invmod(a, n);
        exp = -exp;
    }

    return n_powmod_ui_precomp(a, exp, n, npre);
}
