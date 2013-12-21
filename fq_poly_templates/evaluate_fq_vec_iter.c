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
    Copyright (C) 2012 William Hart
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
_TEMPLATE4(T, poly_evaluate, T, vec_iter)(TEMPLATE(T, struct) * ys,
                                          const TEMPLATE(T, struct) * coeffs, slong len,
                                          const TEMPLATE(T, struct) * xs, slong n,
                                          const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    for (i = 0; i < n; i++)
        _TEMPLATE3(T, poly_evaluate, T)(ys + i, coeffs, len, xs + i, ctx);
}

void
TEMPLATE4(T, poly_evaluate, T, vec_iter)(TEMPLATE(T, struct) * ys,
                                         const TEMPLATE(T, poly_t) poly,
                                         const TEMPLATE(T, struct) * xs, slong n,
                                         const TEMPLATE(T, ctx_t) ctx)
{
    _TEMPLATE4(T, poly_evaluate, T, vec_iter)(ys, poly->coeffs, poly->length,
                                              xs, n, ctx);
}


#endif
