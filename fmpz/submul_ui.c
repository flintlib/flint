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
fmpz_submul_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c1 = *g, r;
    if ((x == 0) || (c1 == 0))
        return;                 /* product is zero */

    r = *f;
    if (r == 0)
    {
        fmpz_mul_ui(f, g, x);   /* we are subtracting product from 0 */
        fmpz_neg(f, f);
        return;
    }

    if (!COEFF_IS_MPZ(c1))      /* c1 is small */
    {
        mp_limb_t prod[2];
        ulong uc1 = FLINT_ABS(c1);

        __mpz_struct * mpz_ptr;
        mpz_t temp;

        umul_ppmm(prod[1], prod[0], uc1, x);    /* compute product */

        if (prod[1] == 0)       /* product fits in one limb */
        {
            if (c1 < 0L)
                fmpz_add_ui(f, f, prod[0]);
            else
                fmpz_sub_ui(f, f, prod[0]);
            return;
        }
        else if ((prod[1] == 1) && (!COEFF_IS_MPZ(r)) && ((r ^ c1) >= 0L))
        {
            /*
               only chance at cancellation is if product is one bit past 
               a limb and f is small and same sign as this product
             */
            ulong ur = FLINT_ABS(r);
            if (ur > prod[0])   /* cancellation will occur */
            {
                fmpz_set_ui(f, prod[0] - ur);
                if (r > 0L)
                    fmpz_neg(f, f);
                return;
            }
        }

        /*
           in all remaining cases res is either big already, 
           or will be big in the end
         */
        mpz_ptr = _fmpz_promote_val(f);
        /* set up a temporary, cheap mpz_t to contain prod */
        temp->_mp_d = prod;
        temp->_mp_size = (c1 < 0L ? -2 : 2);
        mpz_sub(mpz_ptr, mpz_ptr, temp);
        _fmpz_demote_val(f);    /* cancellation may have occurred */
    }
    else                        /* c1 is large */
    {
        __mpz_struct *mpz_ptr = _fmpz_promote_val(f);
        mpz_submul_ui(mpz_ptr, COEFF_TO_PTR(c1), x);
        _fmpz_demote_val(f);    /* cancellation may have occurred */
    }
}
