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

===============================================================================*/
/******************************************************************************

 Copyright (C) 2010 William Hart
 
******************************************************************************/

#ifndef MFPR_POLY_H
#define MPFR_POLY_H

#include <mpir.h>
#include <mpfr.h> 
#include "flint.h"

typedef struct
{
   __mpfr_struct * coeffs;
   ulong length;
   ulong alloc;
   mp_bitcnt_t prec;
} mpfr_poly_struct;

// fmpz_poly_t allows reference-like semantics for fmpz_poly_struct
typedef mpfr_poly_struct mpfr_poly_t[1];

extern gmp_randstate_t mpfr_poly_randstate;

void mpfr_poly_init(mpfr_poly_t poly, mp_bitcnt_t prec);

void mpfr_poly_init2(mpfr_poly_t poly, const ulong alloc, mp_bitcnt_t prec);

void mpfr_poly_realloc(mpfr_poly_t poly, const ulong alloc);

void mpfr_poly_fit_length(mpfr_poly_t poly, const ulong length);

void mpfr_poly_clear(mpfr_poly_t poly);

static inline
void _mpfr_poly_set_length(mpfr_poly_t poly, ulong length)
{
   poly->length = length;
}

void mpfr_poly_randinit(void);

void mpfr_poly_randclear(void);

void mpfr_poly_randtest(mpfr_poly_t poly, ulong length);

#endif






