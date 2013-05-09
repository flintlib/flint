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
fmpz_set_mpz(fmpz_t f, const mpz_t x)
{
    len_t size = (len_t) x->_mp_size;

    if (size == 0L)             /* x is zero */
    {
        fmpz_zero(f);
    }
    else if (size == 1L)        /* x is positive and 1 limb */
    {
        fmpz_set_ui(f, mpz_get_ui(x));
    }
    else if (size == -1L)       /* x is negative and 1 limb */
    {
        ulong uval = mpz_get_ui(x);
        if (uval <= COEFF_MAX)  /* x is small */
        {
            _fmpz_demote(f);
            *f = -uval;
        }
        else                    /* x is large but one limb */
        {
            __mpz_struct *mpz_ptr = _fmpz_promote(f);
            mpz_set_ui(mpz_ptr, uval);
            mpz_neg(mpz_ptr, mpz_ptr);
        }
    }
    else                        /* x is more than one limb */
    {
        __mpz_struct *mpz_ptr = _fmpz_promote(f);
        mpz_set(mpz_ptr, x);
    }
}
