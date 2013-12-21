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

    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
_TEMPLATE3(T, poly_evaluate, T) (TEMPLATE(T, t) rop,
                                 const TEMPLATE(T, struct) * op, slong len,
                                 const TEMPLATE(T, t) a,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    if (len == 0)
    {
        TEMPLATE(T, zero) (rop, ctx);
    }
    else if (len == 1 || TEMPLATE(T, is_zero) (a, ctx))
    {
        TEMPLATE(T, set) (rop, op + 0, ctx);
    }
    else
    {
        slong i = len - 1;
        TEMPLATE(T, t) t;

        TEMPLATE(T, init) (t, ctx);
        TEMPLATE(T, set) (rop, op + i, ctx);
        for (i = len - 2; i >= 0; i--)
        {
            TEMPLATE(T, mul) (t, rop, a, ctx);
            TEMPLATE(T, add) (rop, op + i, t, ctx);
        }
        TEMPLATE(T, clear) (t, ctx);
    }
}

void
TEMPLATE3(T, poly_evaluate, T) (TEMPLATE(T, t) rop,
                                const TEMPLATE(T, poly_t) f,
                                const TEMPLATE(T, t) a,
                                const TEMPLATE(T, ctx_t) ctx)
{
    if (rop == a)
    {
        TEMPLATE(T, t) t;
        TEMPLATE(T, init) (t, ctx);
        _TEMPLATE3(T, poly_evaluate, T) (t, f->coeffs, f->length, a, ctx);
        TEMPLATE(T, swap) (rop, t, ctx);
        TEMPLATE(T, clear) (t, ctx);
    }
    else
    {
        _TEMPLATE3(T, poly_evaluate, T) (rop, f->coeffs, f->length, a, ctx);
    }
}


#endif
