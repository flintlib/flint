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
 Copyright (C) 2010 Daniel Woodhouse
 
******************************************************************************/

#ifndef NMOD_MPOLY_H
#define NMOD_MPOLY_H

#include <stdio.h>
#include <mpir.h>
#include "nmod_vec.h"

typedef struct
{
   mp_ptr coeffs;
   mp_ptr exps;
   long alloc;
   long length;
   ulong vars; // number of variables
   mp_bitcnt_t ebits; // width of exponent bitfield (per variable)
   nmod_t mod;
} nmod_mpoly_struct;

typedef nmod_mpoly_struct nmod_mpoly_t[1];

typedef struct nmod_mpoly_entry_t nmod_mpoly_entry_t;

struct nmod_mpoly_entry_t
{
   nmod_mpoly_entry_t * next; // linked list of entries
   unsigned int i1; // index of coefficient in first poly 
   unsigned int i2; // index of coefficient in second poly
};

typedef struct
{
   mp_limb_t exp;
   nmod_mpoly_entry_t * entry; 
} nmod_mpoly_heap_t;


void nmod_mpoly_init(nmod_mpoly_t poly, mp_limb_t n, long vars, ulong ebits);

void nmod_mpoly_init_preinv(nmod_mpoly_t poly, mp_limb_t n, mp_limb_t ninv, long vars, ulong ebits);

void nmod_mpoly_init2(nmod_mpoly_t poly,mp_limb_t n, long alloc, long vars, ulong ebits);

void nmod_mpoly_init2_preinv(nmod_mpoly_t poly, mp_limb_t n, mp_limb_t ninv, long alloc, long vars, ulong ebits);

void nmod_mpoly_realloc(nmod_mpoly_t poly, long alloc);

void nmod_mpoly_clear(nmod_mpoly_t poly);

void nmod_mpoly_fit_length(nmod_mpoly_t poly, long alloc);

static inline
mp_bitcnt_t nmod_mpoly_max_bits(nmod_mpoly_t poly)
{
   return _nmod_vec_max_bits(poly->coeffs, poly->length);
}


static inline
void _nmod_mpoly_normalise(nmod_mpoly_t poly)
{
   while (poly->length && (poly->coeffs[poly->length - 1] == 0L))
      poly->length--;
}



void nmod_mpoly_randtest(nmod_mpoly_t poly, long length);

void nmod_mpoly_reheapify(nmod_mpoly_heap_t * heap, ulong * n);

void nmod_mpoly_heap_insert(nmod_mpoly_heap_t * heap, ulong * n, nmod_mpoly_entry_t * entry, mp_limb_t exp);

void nmod_mpoly_mul_heap(nmod_mpoly_t res, nmod_mpoly_t poly1, nmod_mpoly_t poly2);


#endif






