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
#include "fmpz.h"

ulong
z_gcd(len_t a, len_t b)
{
    ulong ua = FLINT_ABS(a);
    ulong ub = FLINT_ABS(b);

    return n_gcd_full(ua, ub);
}

void
fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(g))
    {
        fmpz_abs(f, h);
        return;
    }

    if (fmpz_is_zero(h))
    {
        fmpz_abs(f, g);
        return;
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            fmpz_set_si(f, z_gcd(c1, c2));
        }
        else                    /* h is large, but g is small */
        {
            fmpz c2d = fmpz_fdiv_ui(h, FLINT_ABS(c1));
            fmpz_set_si(f, z_gcd(c1, c2d));
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(c2))  /* h is small, but g is large */
        {
            fmpz c1d = fmpz_fdiv_ui(g, FLINT_ABS(c2));
            fmpz_set_si(f, z_gcd(c2, c1d));
        }
        else                    /* g and h are both large */
        {
            __mpz_struct *mpz_ptr = _fmpz_promote(f);   /* aliasing fine as g, h already large */

            mpz_gcd(mpz_ptr, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);    /* gcd may be small */
        }
    }
}
