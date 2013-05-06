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
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_sqrtmod(mp_limb_t a, mp_limb_t p) 
{
    long i, r, m;
    mp_limb_t p1, k, b, g, bpow, gpow, res;
    mp_limb_t pinv;

    if (a <= 1)
    {
        return a;
    }

    pinv = n_preinvert_limb(p);

    if (n_jacobi_unsigned(a, p) == -1)
        return 0;

    if ((p & 3UL) == 3)
    {
        return n_powmod2_ui_preinv(a, (p + 1) / 4, p, pinv);
    }

    r = 0;
    p1 = p - 1;

    do {
        p1 >>= 1UL; 
        r++;
    } while ((p1 & 1UL) == 0);

    b = n_powmod2_ui_preinv(a, p1, p, pinv);

    for (k = 2; ; k++)
    {
        if (n_jacobi_unsigned(k, p) == -1) break;
    }

    g = n_powmod2_ui_preinv(k, p1, p, pinv);
    res = n_powmod2_ui_preinv(a, (p1 + 1) / 2, p, pinv);

    while (b != 1)
    {
        bpow = b;
        m = 0;
        do
        {
            bpow = n_mulmod2_preinv(bpow, bpow, p, pinv);
            m++;
        } while (m < r && bpow != 1);
        gpow = g;
        for (i = 1; i < r - m; i++)
        {
            gpow = n_mulmod2_preinv(gpow, gpow, p, pinv);
        }
        res = n_mulmod2_preinv(res, gpow, p, pinv);
        g = n_mulmod2_preinv(gpow, gpow, p, pinv);
        b = n_mulmod2_preinv(b, g, p, pinv);
        r = m;
    }

    return res;
}

