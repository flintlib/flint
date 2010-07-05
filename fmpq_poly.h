/*============================================================================

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

 Copyright (C) 2010 Sebastian Pancratz
 Copyright (C) 2010 William Hart
 
******************************************************************************/

#ifndef FMPQ_POLY_H
#define FMPQ_POLY_H

#include <mpir.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

typedef struct
{
    fmpz * coeffs;
    fmpz_t den;
    ulong alloc;
    ulong length;
} fmpq_poly_struct;

typedef fmpq_poly_struct fmpq_poly_t[1];

void fmpq_poly_init(fmpq_poly_t poly);

void fmpq_poly_init2(fmpq_poly_t poly);

void fmpq_poly_realloc(fmpq_poly_t poly, const ulong alloc);

void fmpq_poly_clear(fmpq_poly_t poly);

void fmpq_poly_fit_length(fmpq_poly_t poly, ulong length);

void _fmpq_poly_normalise(fmpq_poly_t poly);

void fmpq_poly_canonicalise(fmpq_poly_t poly);

void _fmpq_poly_set_length(fmpq_poly_t poly, const ulong length);

void fmpq_poly_set(fmpq_poly_t poly1, const fmpq_poly_t poly2);

int fmpq_poly_equal(const fmpq_poly_t poly1, const fmpq_poly_t poly2);

static inline 
int fmpq_poly_is_zero(const fmpq_poly_t poly)
{
    return poly->length == 0;
}

#endif

