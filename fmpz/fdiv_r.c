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
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void
fmpz_fdiv_r(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        printf("Exception (fmpz_fdiv_r). Division by zero.\n");
        abort();
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            fmpz q = c1 / c2;   /* compute C quotient */
            fmpz r = c1 - c2 * q;   /* compute remainder */

            if ((c2 > 0L && r < 0L) || (c2 < 0L && r > 0L))
                r += c2;

            fmpz_set_si(f, r);
        }
        else                    /* h is large and g is small */
        {
            if (c1 == 0L)
            {
                fmpz_set_si(f, c1);
            }
            else if ((c1 < 0L && fmpz_sgn(h) < 0) || (c1 > 0L && fmpz_sgn(h) > 0))  /* signs are the same */
            {
                fmpz_set_si(f, c1);
            }
            else
            {
                fmpz_add(f, g, h);
            }
        }
    }
    else                        /* g is large */
    {
        __mpz_struct * mpz_ptr = _fmpz_promote(f);

        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            if (c2 > 0)         /* h > 0 */
            {
                mpz_fdiv_r_ui(mpz_ptr, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                mpz_cdiv_r_ui(mpz_ptr, COEFF_TO_PTR(c1), -c2);
            }
        }
        else                    /* both are large */
        {
            mpz_fdiv_r(mpz_ptr, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
        }
        _fmpz_demote_val(f);    /* division by h may result in small value */
    }
}
