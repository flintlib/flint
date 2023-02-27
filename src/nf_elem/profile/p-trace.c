/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "nf.h"
#include "nf_elem.h"

#define BITS 10

typedef struct
{
   slong length;
   int monic;
} info_t;

void random_fmpq_poly(fmpq_poly_t pol, flint_rand_t state, slong length)
{
   fmpz * arr;
   slong i;

   fmpq_poly_fit_length(pol, length);

   arr = fmpq_poly_numref(pol);

   for (i = 0; i < length; i++)
      fmpz_randbits(arr + i, state, BITS);

   fmpz_randbits(fmpq_poly_denref(pol), state, BITS);

   _fmpq_poly_set_length(pol, length);
   _fmpq_poly_normalise(pol);
   fmpq_poly_canonicalise(pol);
}

void random_nf_elem(nf_elem_t a, flint_rand_t state, nf_t nf)
{
   slong len = nf->pol->length - 1;
   slong i;

   random_fmpq_poly(NF_ELEM(a), state, len);
}

void sample(void * arg, ulong count)
{
   info_t * info = (info_t *) arg;
   slong length = info->length, i, j;
   int monic = info->monic;
   int scale;
   
   scale = 1000;
   if (length >= 50) scale = 100;
   if (length >= 500) scale = 40;
   
   flint_rand_t state;
   flint_randinit(state);

   fmpq_poly_t pol;
   nf_t nf;
   nf_elem_t a;
   fmpq_t norm;
        
   fmpq_poly_init(pol);
   fmpq_init(norm);
     
   for (i = 0; i < count; i++)
   {
      random_fmpq_poly(pol, state, length);
      if (monic)
      {
         fmpz_one(fmpq_poly_denref(pol));
         fmpq_poly_set_coeff_ui(pol, length - 1, 1);
      }
	
      nf_init(nf, pol);
       
      nf_elem_init(a, nf);
        
      random_nf_elem(a, state, nf);
      if (monic)
         fmpz_one(fmpq_poly_denref(NF_ELEM(a)));
	
      prof_start();
      for (j = 0; j < scale; j++)
      {
         nf_elem_trace(norm, a, nf);
      }
      prof_stop();
   }
  
   fmpq_clear(norm);

   nf_elem_clear(a, nf);
        
   nf_clear(nf);

   fmpq_poly_clear(pol);

   flint_randclear(state);
}

int main(void)
{
   double min, max;
   info_t info;
   slong k, scale;

   printf("Number field element trace\n");
   flint_printf("bits = %ld\n", BITS);

   for (k = 4; k <= 1000; k = (slong) ceil(1.1*k))
   {
      info.length = k;
      info.monic = 0;

      scale = 1000;
      if (k >= 50) scale = 100;
      if (k >= 500) scale = 40;
      
      prof_repeat(&min, &max, sample, (void *) &info);
      
      flint_printf("generic: length %wd, min %.3e us, max %.3e us\n", 
           info.length,
		   (min/scale),
           (max/scale)
	     );

      info.monic = 1;
     
      prof_repeat(&min, &max, sample, (void *) &info);
         
      flint_printf("monic  : length %wd, min %.3e us, max %.3e us\n", 
           info.length,
		   (min/scale),
           (max/scale)
	     );
   }

   return 0;
}
