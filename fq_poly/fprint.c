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
#include <stdio.h>
#include <gmp.h>
#include "fmpz.h"
#include "fq_poly.h"

int
_fq_poly_fprint(FILE * file, const fq_struct * poly, slong len,
                const fq_ctx_t ctx)
{
    int r;
    slong i;

    r = fprintf(file, "%ld ", len);
    if (r <= 0)
        return r;

    if (len == 0)
        return r;

    for (i = 0; (r > 0) && (i < len); i++)
    {
        r = fprintf(file, " ");
        if (r <= 0)
            return r;
        r = fq_fprint(file, poly + i, ctx);
        if (r <= 0)
            return r;
    }

    return r;
}

int
fq_poly_fprint(FILE * file, const fq_poly_t poly, const fq_ctx_t ctx)
{
    return _fq_poly_fprint(file, poly->coeffs, poly->length, ctx);
}
