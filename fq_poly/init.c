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

    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "fq_poly.h"

fq_struct * _fq_poly_init(long len)
{
    long i;
    fq_struct *v;

    v = flint_malloc(len * sizeof (fq_struct));

    for (i = 0; i < len; i++)
        fq_init(v + i);

    return v;
}

void fq_poly_init(fq_poly_t poly)
{
    poly->coeffs = NULL;
    poly->alloc  = 0;
    poly->length = 0;
}

void fq_poly_init2(fq_poly_t poly, long alloc)
{
    poly->coeffs = (alloc) ? _fq_poly_init(alloc) : NULL;
    poly->alloc  = alloc;
    poly->length = 0;
}
