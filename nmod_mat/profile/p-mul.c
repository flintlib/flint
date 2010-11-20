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

    Copyright 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

typedef struct
{
   ulong dim;
   mp_limb_t modulus;
} mat_mul_t;

void sample(void * arg, ulong count)
{
   mat_mul_t * params = (mat_mul_t *) arg;
   mp_limb_t n = params->modulus;
   ulong i, dim = params->dim;
   
   nmod_mat_t A, B, C;

   nmod_mat_init(A, dim, dim, n);
   nmod_mat_init(B, dim, dim, n);
   nmod_mat_init(C, dim, dim, n);

   prof_start();
      
   for (i = 0; i < count; i++)
   {
      nmod_mat_mul(C, A, B);
   }

   prof_stop();

   if (C->entries[0] == 5893479483L) abort();

   nmod_mat_clear(A);
   nmod_mat_clear(B);
   nmod_mat_clear(C);
}

int main(void)
{
   double min, max;
   mat_mul_t params;
   long dim;

   printf("nmod_mat_mul:\n");
   
   params.modulus = 40000;
	  
   for (dim = 2; dim <= 512; dim = (long) ((double) dim * 1.1) + 1)
   {
      params.dim = dim;
	  prof_repeat(&min, &max, sample, &params);
      printf("dim = %ld, time is %.2f us\n", 
		  dim, min);
   }
   
   return 0;
}
