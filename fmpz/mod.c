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

void
fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            long r;
            if (c2 < 0L)
                c2 = -c2;
            if (c1 < 0L)
            {
                r = c2 - (-c1 % c2);    /* C doesn't correctly handle negative mods */
                if (r == c2)
                    r = 0;
            }
            else
                r = c1 % c2;

            fmpz_set_si(f, r);
        }
        else                    /* h is large and g is small */
        {
            if (c1 < 0L)
            {
                fmpz_abs(f, h);
                fmpz_sub_ui(f, f, -c1);
            }
            else
                fmpz_set_ui(f, c1);
        }
    }
    else                        /* g is large */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            if (c2 < 0L)
                fmpz_set_si(f, mpz_fdiv_ui(COEFF_TO_PTR(c1), -c2));
            else
                fmpz_set_ui(f, mpz_fdiv_ui(COEFF_TO_PTR(c1), c2));
        }
        else                    /* both are large */
        {
            __mpz_struct *mpz_ptr = _fmpz_promote(f);
            mpz_mod(mpz_ptr, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);    /* reduction mod h may result in small value */
        }
    }
}
