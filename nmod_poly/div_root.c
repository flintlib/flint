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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

mp_limb_t
_nmod_poly_div_root(mp_ptr Q, mp_srcptr A, len_t len, mp_limb_t c, nmod_t mod)
{
    mp_limb_t r, t;
    len_t i;

    if (len < 2)
        return 0UL;

    t = A[len - 2];
    r = Q[len - 2] = A[len - 1];

    for (i = len - 2; i > 0; i--)
    {
        r = nmod_add(nmod_mul(r, c, mod), t, mod);
        t = A[i-1];
        Q[i-1] = r;
    }

    r = nmod_add(nmod_mul(r, c, mod), t, mod);
    return r;
}

mp_limb_t
nmod_poly_div_root(nmod_poly_t Q, 
                 const nmod_poly_t A, mp_limb_t c)
{
    mp_limb_t rem;

    len_t len = A->length;

    if (len == 0)
    {
        nmod_poly_zero(Q);
        return 0UL;
    }

    if (len == 1)
    {
        rem = A->coeffs[0];
        nmod_poly_zero(Q);
        return rem;
    }

    if (c == 0)
    {
        rem = A->coeffs[0];
        nmod_poly_shift_right(Q, A, 1);
        return rem;
    }

    nmod_poly_fit_length(Q, len - 1);
    rem = _nmod_poly_div_root(Q->coeffs, A->coeffs, len, c, Q->mod);
    Q->length = len - 1;
    return rem;
}
