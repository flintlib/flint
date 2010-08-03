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

void nmod_mpoly_init(nmod_mpoly_t poly, mp_limb_t n, long vars, ulong ebits);

void nmod_mpoly_init_preinv(nmod_mpoly_t poly, 
							 mp_limb_t n, mp_limb_t ninv, long vars, ulong ebits);

void nmod_mpoly_init2(nmod_mpoly_t poly, 
					             mp_limb_t n, long alloc, long vars, ulong ebits);

void nmod_mpoly_init2_preinv(nmod_mpoly_t poly, 
			    mp_limb_t n, mp_limb_t ninv, long alloc, long vars, ulong ebits);

void nmod_mpoly_realloc(nmod_mpoly_t poly, long alloc);

void nmod_mpoly_clear(nmod_mpoly_t poly);

void nmod_mpoly_fit_length(nmod_mpoly_t poly, long alloc);

#endif






