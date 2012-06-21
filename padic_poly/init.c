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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "padic_poly.h"

void padic_poly_init(padic_poly_t poly)
{
    poly->coeffs = NULL;
    poly->alloc  = 0;
    poly->length = 0;
    poly->val    = 0;
}

void padic_poly_init2(padic_poly_t poly, long alloc)
{
    poly->coeffs = alloc ? (fmpz *) flint_calloc(alloc, sizeof(fmpz)) : NULL;
    poly->alloc  = alloc;
    poly->length = 0;
    poly->val    = 0;
}

