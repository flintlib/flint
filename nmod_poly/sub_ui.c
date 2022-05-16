/*
    Copyright (C) 2015 Tommy Hofmann
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"


void nmod_poly_sub_ui(nmod_poly_t res, const nmod_poly_t poly, ulong c)
{
    if (c >= poly->mod.n)
        NMOD_RED(c, c, poly->mod);

    if (poly->length == 0)
    {
        if (c == 0)
            nmod_poly_zero(res);
        else
        {
            nmod_poly_fit_length(res, 1);
            nmod_poly_set_coeff_ui(res, 0, poly->mod.n - c);
            _nmod_poly_set_length(res, 1);
        }
    }
    else
    {
        nmod_poly_set(res, poly);

        nmod_poly_set_coeff_ui(res, 0, nmod_sub(res->coeffs[0], c, poly->mod));

        _nmod_poly_normalise(res);

   }
}

