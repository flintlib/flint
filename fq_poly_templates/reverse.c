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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_reverse) (TEMPLATE(T, struct) * res,
                            const TEMPLATE(T, struct) * poly, slong len,
                            slong n, const TEMPLATE(T, ctx_t) ctx)
{
    if (res == poly)
    {
        slong i;

        for (i = 0; i < n / 2; i++)
        {
            TEMPLATE(T, struct) t = res[i];
            res[i] = res[n - 1 - i];
            res[n - 1 - i] = t;
        }

        for (i = 0; i < n - len; i++)
            TEMPLATE(T, zero) (res + i, ctx);
    }
    else
    {
        slong i;

        for (i = 0; i < n - len; i++)
            TEMPLATE(T, zero) (res + i, ctx);

        for (i = 0; i < len; i++)
            TEMPLATE(T, set) (res + (n - len) + i, poly + (len - 1) - i, ctx);
    }
}

void
TEMPLATE(T, poly_reverse) (TEMPLATE(T, poly_t) res,
                           const TEMPLATE(T, poly_t) poly, slong n,
                           const TEMPLATE(T, ctx_t) ctx)
{
    slong len = FLINT_MIN(n, poly->length);
    if (len == 0)
    {
        TEMPLATE(T, poly_zero) (res, ctx);
        return;
    }

    TEMPLATE(T, poly_fit_length) (res, n, ctx);

    _TEMPLATE(T, poly_reverse) (res->coeffs, poly->coeffs, len, n, ctx);

    _TEMPLATE(T, poly_set_length) (res, n, ctx);
    _TEMPLATE(T, poly_normalise) (res, ctx);
}


#endif
