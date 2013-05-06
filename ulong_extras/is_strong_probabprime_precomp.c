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

    Copyright (C) 2008, Peter Shrimpton
    Copyright (C) 2009 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int
n_is_strong_probabprime_precomp(mp_limb_t n, double npre, mp_limb_t a,
                                mp_limb_t d)
{
    mp_limb_t t = d;
    mp_limb_t y;

    y = n_powmod_ui_precomp(a, t, n, npre);

    if (y == 1UL)
        return 1;
    t <<= 1;

    while ((t != n - 1) && (y != n - 1))
    {
        y = n_mulmod_precomp(y, y, n, npre);
        t <<= 1;
    }

    return (y == n - 1);
}
