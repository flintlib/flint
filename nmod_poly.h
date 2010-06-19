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

#include <stdio.h>
#include <mpir.h>
#include "nmod_vec.h"

typedef struct
{
   mp_ptr coeffs;
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

static inline
int nmod_poly_equal(nmod_poly_t a, nmod_poly_t b)
{
   if (a->length != b->length)
	  return 0;

   if (a != b)
      if (!_nmod_vec_equal(a->coeffs, b->coeffs, a->length))
	     return 0;

   return 1;
}

static inline
void nmod_poly_print(nmod_poly_t a)
{
   ulong i;
   
   if (a->length == 0)
   {
	  printf("0");
	  return;
   }
   
   printf("%lu  ", a->length);

   for (i = 0; i < a->length - 1; i++)
      printf("%lu ", a->coeffs[i]);
   printf("%lu", a->coeffs[i]);
}

static inline
void nmod_poly_zero(nmod_poly_t res)
{
   res->length = 0;
}

void nmod_poly_neg(nmod_poly_t res, const nmod_poly_t poly1);

void _nmod_poly_add(mp_ptr res, mp_srcptr poly1, ulong len1, 
					           mp_srcptr poly2, ulong len2, nmod_t mod);

void nmod_poly_add(nmod_poly_t res, const nmod_poly_t poly1, 
				                               const nmod_poly_t poly2);

void _nmod_poly_sub(mp_ptr res, mp_srcptr poly1, ulong len1, 
					           mp_srcptr poly2, ulong len2, nmod_t mod);

void nmod_poly_sub(nmod_poly_t res, const nmod_poly_t poly1, 
				                               const nmod_poly_t poly2);

void nmod_poly_scalar_mul(nmod_poly_t res, 
						          const nmod_poly_t poly1, mp_limb_t c);

#endif






