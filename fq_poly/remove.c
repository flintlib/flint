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

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

ulong
fq_poly_remove(fq_poly_t f, const fq_poly_t g, const fq_ctx_t ctx)
{
    fq_poly_t q, r;
    ulong i = 0;

    fq_poly_init(q, ctx);
    fq_poly_init(r, ctx);

    while (1)
    {
        if (f->length < g->length)
            break;
        fq_poly_divrem(q, r, f, g, ctx);
        if (r->length == 0)
            fq_poly_swap(q, f, ctx);
        else
            break;
        i++;
    }

    fq_poly_clear(q, ctx);
    fq_poly_clear(r, ctx);

    return i;
}
