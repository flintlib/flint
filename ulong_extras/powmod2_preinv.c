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

    Copyright (C) 2009, 2013, 2016 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

ulong
n_powmod2_preinv(ulong a, slong exp, ulong n, ulong ninv)
{
    ulong norm;

    FLINT_ASSERT(n != 0);

    if (a >= n)
        a = n_mod2_preinv(a, n, ninv);

    if (exp < 0)
    {
        ulong g = n_gcdinv(&a, a, n);

        if (g != 1)
            flint_throw(FLINT_IMPINV, "Cannot invert modulo %wd*%wd\n", g,
                        n / g);

        exp = -exp;
    }

    count_leading_zeros(norm, n);

    return n_powmod_ui_preinv(a << norm, exp, n << norm, ninv, norm) >> norm;
}
