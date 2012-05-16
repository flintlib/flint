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

   Copyright (C) 2012 Thomas M. DuBuisson

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_or(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
        fmpz c1,c2;
        fmpz_t tmp;
        fmpz_init(tmp);
        c1 = *g;
        c2 = *h;
        if(!COEFF_IS_MPZ(c1))
        {
            if(!COEFF_IS_MPZ(c2)) /* both inputs are small */
            {
                fmpz_set_si(tmp, c1 | c2);
            } else /* g is small, h is large */
            {
                __mpz_struct * mpz3 = _fmpz_promote(tmp);
                __mpz_struct * mpz2 = COEFF_TO_PTR(c2);
                mpz_t g2;
                mpz_init_set_si(g2,c1);
                mpz_ior(mpz3, mpz2, g2);
                mpz_clear(g2);
            }
        } else
        {
            if(!COEFF_IS_MPZ(c2)) /* g is large, h is small */
            {
                __mpz_struct *mpz3 = _fmpz_promote(tmp);
                __mpz_struct *mpz1 = COEFF_TO_PTR(c1);
                mpz_t h2;
                mpz_init_set_si(h2,c2);
                mpz_ior(mpz3, mpz1, h2);
                mpz_clear(h2);
            } else /* g and h are large */
            {
                __mpz_struct * mpz3 = _fmpz_promote(tmp);
                __mpz_struct * mpz1 = COEFF_TO_PTR(c1);
                __mpz_struct * mpz2 = COEFF_TO_PTR(c2);
                mpz_ior(mpz3, mpz1, mpz2);
            }
        }
        fmpz_set(f,tmp);
        fmpz_clear(tmp);
}

