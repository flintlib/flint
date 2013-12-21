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
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_inflate) (TEMPLATE(T, poly_t) result,
                           const TEMPLATE(T, poly_t) input, ulong inflation,
                           const TEMPLATE(T, ctx_t) ctx)
{
    if (input->length <= 1 || inflation == 1)
    {
        TEMPLATE(T, poly_set) (result, input, ctx);
    }
    else if (inflation == 0)
    {
        TEMPLATE(T, t) v;
        TEMPLATE(T, init) (v, ctx);
        TEMPLATE(T, one) (v, ctx);
        TEMPLATE(T, TEMPLATE(poly_evaluate, T)) (v, input, v, ctx);
        TEMPLATE(T, poly_zero) (result, ctx);
        TEMPLATE(T, poly_set_coeff) (result, 0, v, ctx);
        TEMPLATE(T, clear) (v, ctx);
    }
    else
    {
        slong i, j, res_length = (input->length - 1) * inflation + 1;

        TEMPLATE(T, poly_fit_length) (result, res_length, ctx);

        for (i = input->length - 1; i > 0; i--)
        {
            TEMPLATE(T, set) (result->coeffs + (i * inflation),
                              input->coeffs + i, ctx);
            for (j = i * inflation - 1; j > (i - 1) * inflation; j--)
                TEMPLATE(T, zero) (result->coeffs + j, ctx);
        }
        TEMPLATE(T, set) (result->coeffs, input->coeffs, ctx);
        result->length = res_length;
    }
}


#endif
