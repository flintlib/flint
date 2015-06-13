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

    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"

int
fmpz_get_mpn(mp_ptr *n, fmpz_t n_in)
{
    mp_limb_t n_size;
    mp_ptr temp;

    n_size = fmpz_size(n_in);
    *n = flint_malloc(n_size * sizeof(mp_limb_t));

    if (n_size <= 1)
    {
        (*n)[0] = fmpz_get_ui(n_in);
        return 1;
    }
    else
    {
        temp = COEFF_TO_PTR(*n_in)->_mp_d;
        mpn_copyi(*n, temp, n_size);
        return n_size;
    }
}
