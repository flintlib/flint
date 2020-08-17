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

int fmpz_cmpabs(const fmpz_t f, const fmpz_t g)
{
    if (f == g) return 0;  /* aliased inputs */

    if (!COEFF_IS_MPZ(*f))
    {
        if (!COEFF_IS_MPZ(*g)) 
        {
            mp_limb_t uf = FLINT_ABS(*f);
            mp_limb_t ug = FLINT_ABS(*g);
            
            return (uf < ug ? -1 : (uf > ug));
        }
        else return -1;
    }
    else 
    {
        if (!COEFF_IS_MPZ(*g)) return 1;  /* f is large, so if g isn't... */
        else return mpz_cmpabs(COEFF_TO_PTR(*f), COEFF_TO_PTR(*g)); 
    }
}
