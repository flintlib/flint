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

#include "fq_poly.h"

void
_fq_poly_reverse(fq_struct * res, const fq_struct * poly, slong len, slong n,
                 const fq_ctx_t ctx)
{
    if (res == poly)
    {
        slong i;

        for (i = 0; i < n / 2; i++)
        {
            fq_struct t = res[i];
            res[i] = res[n - 1 - i];
            res[n - 1 - i] = t;
        }

        for (i = 0; i < n - len; i++)
            fq_zero(res + i, ctx);
    }
    else
    {
        slong i;

        for (i = 0; i < n - len; i++)
            fq_zero(res + i, ctx);

        for (i = 0; i < len; i++)
            fq_set(res + (n - len) + i, poly + (len - 1) - i, ctx);
    }
}

void
fq_poly_reverse(fq_poly_t res, const fq_poly_t poly, slong n,
                const fq_ctx_t ctx)
{
    slong len = FLINT_MIN(n, poly->length);
    if (len == 0)
    {
        fq_poly_zero(res, ctx);
        return;
    }

    fq_poly_fit_length(res, n, ctx);

    _fq_poly_reverse(res->coeffs, poly->coeffs, len, n, ctx);

    _fq_poly_set_length(res, n, ctx);
    _fq_poly_normalise(res, ctx);
}
