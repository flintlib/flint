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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

mp_limb_t
_nmod_poly_evaluate_nmod(mp_srcptr poly, long len, mp_limb_t c, nmod_t mod)
{
    long m;
    mp_limb_t val;

    if (len == 0)
        return 0;

    if (len == 1 || c == 0)
        return poly[0];

    m = len - 1;
    
    val = poly[m];
    m--;

    for ( ; m >= 0; m--)
    {
        val = n_mulmod2_preinv(val, c, mod.n, mod.ninv);
        val = n_addmod(val, poly[m], mod.n);
    }

    return val;
}

mp_limb_t
nmod_poly_evaluate_nmod(const nmod_poly_t poly, mp_limb_t c)
{
    return _nmod_poly_evaluate_nmod(poly->coeffs, poly->length, c, poly->mod);
}

