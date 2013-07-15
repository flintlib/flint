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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        printf("Exception (fmpz_divexact). Division by zero.\n");
        abort();
    }

    if (!COEFF_IS_MPZ(c1))  /* g is small, h must be also or division isn't exact */
    {
        fmpz_set_si(f, c1 / c2);
    }
    else  /* g is large */
    {
        __mpz_struct * mpz_ptr = _fmpz_promote(f);

        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            if (c2 > 0)  /* h > 0 */
            {
                mpz_divexact_ui(mpz_ptr, COEFF_TO_PTR(c1), c2);
                _fmpz_demote_val(f);  /* division by h may result in small value */
            }
            else
            {
                mpz_divexact_ui(mpz_ptr, COEFF_TO_PTR(c1), -c2);
                _fmpz_demote_val(f);  /* division by h may result in small value */

                fmpz_neg(f, f);
            }
        }
        else  /* both are large */
        {
            mpz_divexact(mpz_ptr, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);  /* division by h may result in small value */
        }
    }
}
