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
TEMPLATE(T, poly_deflate) (TEMPLATE(T, poly_t) result,
                           const TEMPLATE(T, poly_t) input, ulong deflation,
                           const TEMPLATE(T, ctx_t) ctx)
{
    slong res_length, i;

    if (deflation == 0)
    {
        TEMPLATE_PRINTF("Exception (%s_poly_deflate). Division by zero.\n", T);
        abort();
    }

    if (input->length <= 1 || deflation == 1)
    {
        TEMPLATE(T, poly_set) (result, input, ctx);
        return;
    }

    res_length = (input->length - 1) / deflation + 1;
    TEMPLATE(T, poly_fit_length) (result, res_length, ctx);
    for (i = 0; i < res_length; i++)
        TEMPLATE(T, set) (result->coeffs + i, input->coeffs + (i * deflation),
                          ctx);

    result->length = res_length;
}


#endif
