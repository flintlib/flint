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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "nmod_poly.h"

void
_nmod_poly_taylor_shift_horner(mp_ptr poly, mp_limb_t c, len_t n, nmod_t mod)
{
    len_t i, j;

    if (c == 1)
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                poly[j] = nmod_add(poly[j], poly[j + 1], mod);
    }
    else if (c == mod.n - 1)
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                poly[j] = nmod_sub(poly[j], poly[j + 1], mod);
    }
    else if (c != 0)
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                NMOD_ADDMUL(poly[j], poly[j + 1], c, mod);
    }
}

void
nmod_poly_taylor_shift_horner(nmod_poly_t g, const nmod_poly_t f, mp_limb_t c)
{
    if (f != g)
        nmod_poly_set(g, f);

    _nmod_poly_taylor_shift_horner(g->coeffs, c, g->length, g->mod);
}
