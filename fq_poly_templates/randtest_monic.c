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

    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_randtest_monic) (TEMPLATE(T, poly_t) f, flint_rand_t state,
                                  slong len, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    TEMPLATE(T, poly_fit_length) (f, len, ctx);
    for (i = 0; i < len - 1; i++)
    {
        TEMPLATE(T, randtest) (f->coeffs + i, state, ctx);
    }
    TEMPLATE(T, one) (f->coeffs + (len - 1), ctx);
    _TEMPLATE(T, poly_set_length) (f, len, ctx);
    _TEMPLATE(T, poly_normalise) (f, ctx);
}


#endif
