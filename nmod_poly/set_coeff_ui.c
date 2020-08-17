/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void nmod_poly_set_coeff_ui(nmod_poly_t poly, slong j, ulong c)
{
    if (c >= poly->mod.n)
        NMOD_RED(c, c, poly->mod);

    nmod_poly_fit_length(poly, j + 1);

    if (j + 1 < poly->length) /* interior */
        poly->coeffs[j] = c;
    else if (j + 1 == poly->length) /* leading coeff */
    {
        if (c != 0)
            poly->coeffs[j] = c;
        else
        {
            poly->length--;
            _nmod_poly_normalise(poly);
        }
    } else /* extend polynomial */
    {
        if (c == 0) return;
        else
        {
            flint_mpn_zero(poly->coeffs + poly->length, j - poly->length);

            poly->coeffs[j] = c;
            poly->length = j + 1;
        }
    }
}
