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

    Copyright (C) 2009, 2010, 2012 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t
n_mulmod_preinv(mp_limb_t a, mp_limb_t b, mp_limb_t n, 
                                          mp_limb_t ninv, ulong norm)
{
    mp_limb_t q0, q1, r, p_hi, p_lo;

    /* multiply */
    umul_ppmm(p_hi, p_lo, a, b);
    
    /* renormalise product */
    if (norm)
    {
       p_lo = (p_lo >> norm) + (p_hi << (FLINT_BITS - norm));
       p_hi = (p_hi >> norm);
    }

    /* reduce mod n */
    {
        umul_ppmm(q1, q0, ninv, p_hi);
        add_ssaaaa(q1, q0, q1, q0, p_hi, p_lo);

        r = (p_lo - (q1 + 1) * n);

        if (r >= q0)
            r += n;

        return (r < n ? r : r - n);
    }
}
