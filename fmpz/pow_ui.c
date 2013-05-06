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
fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz c1;

    if (exp == 0L)
    {
        fmpz_one(f);
        return;
    }

    c1 = *g;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        ulong u1 = FLINT_ABS(c1);
        ulong bits = FLINT_BIT_COUNT(u1);
        if ((bits <= 1) || (exp * bits <= FLINT_BITS - 2))
        {
            fmpz_set_ui(f, n_pow(u1, exp));
        }
        else
        {
            __mpz_struct *mpz_ptr = _fmpz_promote_val(f);

            mpz_set_ui(mpz_ptr, u1);
            mpz_pow_ui(mpz_ptr, mpz_ptr, exp);
            _fmpz_demote_val(f);    /* may actually fit into a small after all */
        }

        if ((c1 < 0L) && (exp & 1)) /* sign is -ve if exp odd and g -ve */
            fmpz_neg(f, f);
    }
    else
    {
        __mpz_struct *mpz_ptr = _fmpz_promote_val(f);
        mpz_pow_ui(mpz_ptr, COEFF_TO_PTR(c1), exp);
        /* no need to demote as it can't get smaller */
    }
}
