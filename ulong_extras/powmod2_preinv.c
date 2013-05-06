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
n_powmod2_ui_preinv(mp_limb_t a, mp_limb_t exp, mp_limb_t n, mp_limb_t ninv)
{
    mp_limb_t x, y;

    if (n == 1UL) return 0UL;

    x = 1UL;
    y = a;
    while (exp)
    {
        if (exp & 1) x = n_mulmod2_preinv(x, y, n, ninv);
        exp >>= 1;
        if (exp) y = n_mulmod2_preinv(y, y, n, ninv);
    }

    return x;
}

mp_limb_t
n_powmod2_preinv(mp_limb_t a, mp_limb_signed_t exp, mp_limb_t n, mp_limb_t ninv)
{
    if (exp < 0L)
    {
        a = n_invmod(a, n);
        exp = -exp;
    }

    return n_powmod2_ui_preinv(a, exp, n, ninv);
}

