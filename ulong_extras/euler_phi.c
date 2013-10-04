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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_euler_phi(mp_limb_t n)
{
    int i;
    mp_limb_t phi;
    n_factor_t fac;

    if (n < 2)
        return n;

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    phi = UWORD(1);
    for (i = 0; i < fac.num; i++)
        phi *= (fac.p[i]-1) * n_pow(fac.p[i], fac.exp[i]-1);

    return phi;
}
