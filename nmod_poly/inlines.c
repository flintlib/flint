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

    Copyright (C) 2015 Tommy Hofmann
    Copyright (C) 2015 William Hart

******************************************************************************/

#define NMOD_POLY_INLINES_C

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

void nmod_poly_add_ui(nmod_poly_t res, const nmod_poly_t poly, ulong c)
{
   
    if (poly->length == 0)
    {
        if (c == 0)
            nmod_poly_zero(res);
        else
        {
            nmod_poly_fit_length(res, 1);
            nmod_poly_set_coeff_ui(res, 0, c);
            _nmod_poly_set_length(res, 1);
        }
    }
    else
    {
        if (c >= poly->mod.n)
            NMOD_RED(c, c, poly->mod);

        nmod_poly_set(res, poly);

        nmod_poly_set_coeff_ui(res, 0, nmod_add(res->coeffs[0], c, poly->mod));

        _nmod_poly_normalise(res);
   }
}

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
            nmod_poly_set_coeff_ui(res, 0, c);
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

