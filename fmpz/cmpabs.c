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
