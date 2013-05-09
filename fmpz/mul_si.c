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
fmpz_mul_si(fmpz_t f, const fmpz_t g, len_t x)
{
    fmpz c2 = *g;

    if (x == 0)
    {
        fmpz_zero(f);
        return;
    }
    else if (!COEFF_IS_MPZ(c2)) /* c2 is small */
    {
        mp_limb_t prod[2];
        mp_limb_t uc2 = FLINT_ABS(c2);
        mp_limb_t ux = FLINT_ABS(x);

        /* unsigned limb by limb multiply (assembly for most CPU's) */
        umul_ppmm(prod[1], prod[0], uc2, ux);
        if ((c2 ^ x) < 0L)
            fmpz_neg_uiui(f, prod[1], prod[0]);
        else
            fmpz_set_uiui(f, prod[1], prod[0]);
    }
    else                        /* c2 is large */
    {
        __mpz_struct *mpz_ptr = _fmpz_promote(f);   /* ok without val as if aliased both are large */
        mpz_mul_si(mpz_ptr, COEFF_TO_PTR(c2), x);
    }
}
