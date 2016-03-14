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

/* a = (b + c) mod n 

    Not a normal add mod function, assumes n is normalized (higest bit set)
    and b and c are reduced modulo n

*/

void 
fmpz_factor_ecm_addmod(mp_ptr a, mp_ptr b, mp_ptr c, mp_ptr n, mp_limb_t n_size)
{
	mp_limb_t cy;

    cy = mpn_add_n(a, b, c, n_size);

    if (cy || (mpn_cmp(a, n, n_size) > 0))
        mpn_sub_n(a, a, n, n_size);
}
