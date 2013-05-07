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
fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz d = *g;

    if (!COEFF_IS_MPZ(d))       /* g is small */
    {
        ulong dabs = FLINT_ABS(d);
        ulong bits = FLINT_BIT_COUNT(dabs);
        if (bits == 0)
        {
            fmpz_set_si(f, 0);
        }
        else if (bits + exp <= FLINT_BITS - 2)  /* result will fit in small */
        {
            fmpz_set_si(f, d << exp);
        }
        else                    /* result is large */
        {
            __mpz_struct *mpz_ptr = _fmpz_promote(f);   /* g is saved */
            mpz_set_si(mpz_ptr, d);
            mpz_mul_2exp(mpz_ptr, mpz_ptr, exp);
        }
    }
    else                        /* g is large */
    {
        __mpz_struct *mpz_ptr = _fmpz_promote(f);   /* g is already large */
        mpz_mul_2exp(mpz_ptr, COEFF_TO_PTR(d), exp);
    }
}
