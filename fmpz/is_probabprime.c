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

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
fmpz_is_probabprime(const fmpz_t p)
{
    fmpz c = *p;
    
    if (fmpz_sgn(p) <= 0)
       return 0;

    if (!COEFF_IS_MPZ(c))
       return n_is_probabprime(c);
    else
    {
       int ret;
       
       ret = (mpz_probab_prime_p(COEFF_TO_PTR(c), 25) != 0);
       return ret;
    }
}
