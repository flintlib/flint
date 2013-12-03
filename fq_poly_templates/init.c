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
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_init) (TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx)
{
    poly->coeffs = NULL;
    poly->alloc = 0;
    poly->length = 0;
}

void
TEMPLATE(T, poly_init2) (TEMPLATE(T, poly_t) poly, slong alloc,
                         const TEMPLATE(T, ctx_t) ctx)
{
    poly->coeffs = (alloc) ? _TEMPLATE(T, vec_init) (alloc, ctx) : NULL;
    poly->alloc = alloc;
    poly->length = 0;
}


#endif
