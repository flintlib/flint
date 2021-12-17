/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_addmul_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c1, r;
    
    c1 = *g;
    if ((x == 0) || (c1 == 0))  /* product is zero */
        return;

    r = *f;
    if (r == 0)
    {
        fmpz_mul_ui(f, g, x);  /* we are adding product to 0 */
        return;
    }

    if (x == UWORD(1))  /* special case, adding g*1 to f */
    {
        fmpz_add(f, f, g);
        return;
    }

    if (c1 == UWORD(1))  /* special case, adding 1*x to f */
    {
        fmpz_add_ui(f, f, x);
        return;
    }

    if (!COEFF_IS_MPZ(c1))  /* c1 is small */
    {
        mp_limb_t prod[2];
        ulong uc1 = FLINT_ABS(c1);

        umul_ppmm(prod[1], prod[0], uc1, x);  /* compute product */

        if (prod[1] == 0)  /* product fits in one limb */
        {
            if (c1 < WORD(0))
                fmpz_sub_ui(f, f, prod[0]);
            else
                fmpz_add_ui(f, f, prod[0]);
            return;
        }
        else if ((prod[1] == 1) && (!COEFF_IS_MPZ(r)) && ((r ^ c1) < WORD(0)))
        {
            /*
               only chance at cancellation is if product is one bit past 
               a limb and res is small and opposite sign to this product
             */
            ulong ur = FLINT_ABS(r);
            if (ur > prod[0])  /* cancellation will occur */
            {
                fmpz_set_ui(f, prod[0] - ur);
                if (r > WORD(0))
                    fmpz_neg(f, f);
                return;
            }
        }
        
        /*
           in all remaining cases res is either big already, 
           or will be big in the end
         */
        {
            __mpz_struct * mf = _fmpz_promote_val(f);
            mpz_t temp;  /* set up a temporary, cheap mpz_t to contain prod */
            temp->_mp_d = prod;
            temp->_mp_size = (c1 < WORD(0) ? -2 : 2);
            mpz_add(mf, mf, temp);
            _fmpz_demote_val(f);  /* cancellation may have occurred */
        }
    }
    else  /* c1 is large */
    {
        __mpz_struct * mf = _fmpz_promote_val(f);
        flint_mpz_addmul_ui(mf, COEFF_TO_PTR(c1), x);
        _fmpz_demote_val(f);  /* cancellation may have occurred */
    }
}
