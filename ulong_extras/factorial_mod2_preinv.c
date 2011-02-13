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

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"


mp_limb_t n_factorial_mod2_preinv(ulong n, mp_limb_t p, mp_limb_t pinv)
{
    mp_limb_t prod, hi, lo;

    prod = 1UL;
    lo = 1UL;

    /* TODO: speedup for n in the range of sqrt(ULONG_MAX) */
    while (n)
    {
        umul_ppmm(hi, lo, lo, n);

        if (hi)
        {
            lo = n_ll_mod_preinv(hi, lo, p, pinv);
            prod = n_mulmod2_preinv(prod, lo, p, pinv);
            lo = 1UL;
        }

        n--;
    }

    return n_mulmod2_preinv(prod, lo, p, pinv);
}
