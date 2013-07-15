/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License,or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not,write to the Free Software
    Foundation,Inc.,51 Franklin St,Fifth Floor,Boston,MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"


void n_nth_prime_bounds(mp_limb_t *lo, mp_limb_t *hi, ulong n)
{
    int bits, ll;
    double llo, lhi;

    /* Lower and upper bounds for ln(n) */
    bits = FLINT_BIT_COUNT(n);
    llo = (bits-1) * 0.6931471;
    lhi = bits * 0.6931472;

    /* Lower bound for ln(ln(n)) */
    if      (n < 16)        ll = 0;
    else if (n < 1619)      ll = 1;
    else if (n < 528491312) ll = 2;
    else                    ll = 3;

    *lo = (mp_limb_t) (n * (llo + ll - 1));
    *hi = (mp_limb_t) (n * (lhi + (ll+1) - (n >= 15985 ? 0.9427 : 0.0)));
}
