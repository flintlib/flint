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
#include "fmpz.h"

typedef struct
{
   mp_ptr coeffs;
   mp_ptr exps;
   long alloc;
   long length;
   ulong vars; /* number of variables */
   mp_bitcnt_t ebits; /* width of exponent bitfield (per variable) */
   nmod_t mod;
} nmod_mpoly_struct;

typedef nmod_mpoly_struct nmod_mpoly_t[1];

typedef struct nmod_mpoly_entry_t nmod_mpoly_entry_t;

struct nmod_mpoly_entry_t
{
   nmod_mpoly_entry_t * next; /* linked list of entries */
   unsigned int i1; /* index of coefficient in first poly */
   unsigned int i2; /* index of coefficient in second poly */
};

typedef struct
{
   mp_limb_t exp;
   nmod_mpoly_entry_t * entry; 
} nmod_mpoly_heap_t;


void nmod_mpoly_init(nmod_mpoly_t poly, mp_limb_t n, 
                                       long vars, ulong ebits);

void nmod_mpoly_init_preinv(nmod_mpoly_t poly, mp_limb_t n, 
                       mp_limb_t ninv, long vars, ulong ebits);

void nmod_mpoly_init2(nmod_mpoly_t poly, mp_limb_t n, 
                           long alloc, long vars, ulong ebits);

void nmod_mpoly_init2_preinv(nmod_mpoly_t poly, mp_limb_t n, 
           mp_limb_t ninv, long alloc, long vars, ulong ebits);

void nmod_mpoly_realloc(nmod_mpoly_t poly, long alloc);

void nmod_mpoly_clear(nmod_mpoly_t poly);

void nmod_mpoly_fit_length(nmod_mpoly_t poly, long alloc);

static __inline__
mp_bitcnt_t nmod_mpoly_max_bits(nmod_mpoly_t poly)
{
   return _nmod_vec_max_bits(poly->coeffs, poly->length);
}


static __inline__
void _nmod_mpoly_normalise(nmod_mpoly_t poly)
{
   while (poly->length && (poly->coeffs[poly->length - 1] == 0L))
      poly->length--;
}


/*
   Assumes that the width of the exponent bitfield is the same.
*/
static __inline__
int nmod_mpoly_equal(nmod_mpoly_t a, nmod_mpoly_t b)
{
   if (a->length != b->length)
	  return 0;

   if (a != b)
      if (!_nmod_vec_equal(a->coeffs, b->coeffs, a->length))
	     return 0;
      if (!_nmod_vec_equal(a->exps, b->exps, a->length))
	     return 0;

   return 1;
}

static __inline__
void nmod_mpoly_set(nmod_mpoly_t a, nmod_mpoly_t b)
{
   if (a != b)
   {
      nmod_mpoly_fit_length(a, b->length);
	   mpn_copyi(a->coeffs, b->coeffs, b->length);
      mpn_copyi(a->exps, b->exps, b->length);
	   a->length = b->length;
      a->ebits = b->ebits;
      a->mod = b->mod;
   }
}

void nmod_mpoly_randtest(nmod_mpoly_t poly, flint_rand_t state, long length);

void nmod_mpoly_reheapify(nmod_mpoly_heap_t * heap, ulong * n);

void nmod_mpoly_heap_insert(nmod_mpoly_heap_t * heap, 
               ulong * n, nmod_mpoly_entry_t * entry, mp_limb_t exp);

void nmod_mpoly_mul_heap(nmod_mpoly_t res, 
                             nmod_mpoly_t poly1, nmod_mpoly_t poly2);

ulong nmod_mpoly_get_coeff(nmod_mpoly_t poly, ulong exp);

void get_period_sequence(fmpz_t *zeroCoefficients, long *coefficients, 
         ulong *exponents, ulong length, ulong monomial, 
                int pow, ulong *primes, int numOfPrimes, ulong nvars);

void get_period_sequence2(fmpz_t *zeroCoefficients, long *coefficients,
         ulong *exponents, ulong length, ulong monomial, 
                int pow, ulong *primes, int numOfPrimes, ulong nvars);

ulong nmod_mpoly_get_coeff_of_product(nmod_mpoly_t poly1, 
                                       nmod_mpoly_t poly2, ulong mon);

#endif






