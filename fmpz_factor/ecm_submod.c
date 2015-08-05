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
#include "fmpz.h"
#include "mpn_extras.h"

/* x = (a - b) mod n 

    Not a normal sub mod function, assumes n is normalized (higest bit set)
    and a and b are reduced modulo n

*/

void
fmpz_factor_ecm_submod(mp_ptr x, mp_ptr a, mp_ptr b, mp_ptr n, mp_limb_t n_size)
{
    mp_ptr temp;

    TMP_INIT;

    TMP_START;
    temp = TMP_ALLOC(n_size * sizeof(mp_limb_t));

    if (mpn_cmp(a, b, n_size) > 0)
        mpn_sub_n(x, a, b, n_size);
    else
    {
        mpn_sub_n(temp, n, b, n_size);
        mpn_add_n(x, temp, a, n_size);
    }

    TMP_END;
}
