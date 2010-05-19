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

#ifndef NMOD_POLY_H
#define NMOD_POLY_H

#include <mpir.h>
#include "nmod_vec.h"

typedef struct
{
   mp_limb_t * coeffs;
   ulong alloc;
   ulong length;
   nmod_t mod;
} nmod_poly_struct;

typedef nmod_poly_struct nmod_poly_t[1];

void nmod_poly_init(nmod_poly_t poly, mp_limb_t n);

void nmod_poly_init_preinv(nmod_poly_t poly, 
							 mp_limb_t n, mp_limb_t ninv);

void nmod_poly_init2(nmod_poly_t poly, 
					             mp_limb_t n, ulong alloc);

void nmod_poly_init2_preinv(nmod_poly_t poly, 
			    mp_limb_t n, mp_limb_t ninv, ulong alloc);

void nmod_poly_realloc(nmod_poly_t poly, ulong alloc);

void nmod_poly_clear(nmod_poly_t poly);

void nmod_poly_fit_length(nmod_poly_t poly, ulong alloc);

static inline
void _nmod_poly_normalise(nmod_poly_t poly)
{
   while (poly->length && (poly->coeffs[poly->length - 1] == 0L))
      poly->length--;
}

void nmod_poly_randtest(nmod_poly_t poly, ulong length);

#endif






