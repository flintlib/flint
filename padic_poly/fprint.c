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

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "padic_poly.h"

int _padic_poly_fprint(FILE *file, const fmpz *poly, long val, long len, 
                       const padic_ctx_t ctx)
{
    long i, v;
    fmpz_t u;

    if (len == 0)
    {
        fprintf(file, "0");
        return 1;
    }

    fmpz_init(u);

    fprintf(file, "%ld ", len);

    for (i = 0; i < len; i++)
    {
        fprintf(file, " ");

        if (fmpz_is_zero(poly + i))
        {
            fprintf(file, "0");
        }
        else
        {
            v = val + fmpz_remove(u, poly + i, ctx->p);

            _padic_fprint(file, u, v, ctx);
        }
    }

    fmpz_clear(u);

    return 1;
}

int padic_poly_fprint(FILE *file, const padic_poly_t poly, 
                      const padic_ctx_t ctx)
{
    _padic_poly_fprint(file, poly->coeffs, poly->val, poly->length, ctx);

    return 1;
}

