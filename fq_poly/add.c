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

******************************************************************************/

#include "fq_poly.h"

void _fq_poly_add(fq_struct *res, 
                  const fq_struct *poly1, long len1, 
                  const fq_struct *poly2, long len2, 
                  const fq_ctx_t ctx)
{
    const long min  = FLINT_MIN(len1, len2);
    long i;

    for (i = 0; i < min; i++)
        fq_add(res + i, poly1 + i, poly2 + i, ctx);

    if (poly1 != res)
        for (i = min; i < len1; i++)
            fq_set(res + i, poly1 + i);

    if (poly2 != res)
        for (i = min; i < len2; i++)
            fq_set(res + i, poly2 + i);
}

void fq_poly_add(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, 
                 const fq_ctx_t ctx)
{
    const long max  = FLINT_MAX(poly1->length, poly2->length);

    fq_poly_fit_length(res, max);

    _fq_poly_add(res->coeffs, poly1->coeffs, poly1->length, 
                              poly2->coeffs, poly2->length, ctx);

    _fq_poly_set_length(res, max);
    _fq_poly_normalise(res);
}
