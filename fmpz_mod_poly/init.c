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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_init(fmpz_mod_poly_t poly, const fmpz_t p)
{
    poly->coeffs = NULL;
    poly->alloc  = 0;
    poly->length = 0;
    fmpz_init(&(poly->p));
    fmpz_set(&(poly->p), p);
}

void fmpz_mod_poly_init2(fmpz_mod_poly_t poly, const fmpz_t p, len_t alloc)
{
    if (alloc)                  /* allocate space for alloc small coeffs */
        poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
    else
        poly->coeffs = NULL;

    poly->alloc = alloc;
    poly->length = 0;
    fmpz_init(&(poly->p));
    fmpz_set(&(poly->p), p);
}

