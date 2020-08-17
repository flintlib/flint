/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qsieve.h"

mp_limb_t qsieve_poly_init(qs_t qs_inf)
{
   ulong num_primes = qs_inf->num_primes;
   ulong s = qs_inf->s;    /* number of prime factors in A coeff */
   mp_limb_t ** A_inv2B;
   slong i;

   fmpz_init(qs_inf->A);
   fmpz_init(qs_inf->B);
   fmpz_init(qs_inf->upp_bound);
   fmpz_init(qs_inf->low_bound);

   qs_inf->curr_subset = flint_malloc(s * sizeof(mp_limb_t));
   qs_inf->first_subset = flint_malloc(s * sizeof(mp_limb_t));
   qs_inf->B_terms = flint_malloc(s * sizeof(mp_limb_t));
   qs_inf->A_ind = flint_malloc(s * sizeof(mp_limb_t));
   qs_inf->A_divp = flint_malloc(s * sizeof(mp_limb_t));
   qs_inf->B0_terms = flint_malloc(s * sizeof(mp_limb_t));

   qs_inf->A_inv2B = flint_malloc(s * sizeof(mp_limb_t *));

   qs_inf->A_inv = flint_malloc(num_primes * sizeof(mp_limb_t));
   qs_inf->soln1 = flint_malloc(num_primes * sizeof(mp_limb_t));
   qs_inf->soln2 = flint_malloc(num_primes * sizeof(mp_limb_t));

   qs_inf->poly = flint_malloc((qs_inf->num_handles + 1)* sizeof(qs_poly_s));

   for (i = 0; i <= qs_inf->num_handles ; i++)
   {
      fmpz_init(qs_inf->poly[i].B);
      qs_inf->poly[i].posn1 = flint_malloc((num_primes + 16)*sizeof(mp_limb_t));
      qs_inf->poly[i].posn2 = flint_malloc((num_primes + 16)*sizeof(mp_limb_t));
      qs_inf->poly[i].soln1 = flint_malloc((num_primes + 16)*sizeof(mp_limb_t));
      qs_inf->poly[i].soln2 = flint_malloc((num_primes + 16)*sizeof(mp_limb_t));
      qs_inf->poly[i].small = flint_malloc(qs_inf->small_primes*sizeof(mp_limb_t));
      qs_inf->poly[i].factor = flint_malloc(qs_inf->max_factors*sizeof(fac_t));
   }

   A_inv2B = qs_inf->A_inv2B;

   for (i = 0; i < s; i++)
      A_inv2B[i] = flint_malloc(num_primes * sizeof(mp_limb_t));

   for (i = 0; i < s; i++)
   {
       fmpz_init(qs_inf->A_divp[i]);
       fmpz_init(qs_inf->B_terms[i]);
   }

   return 0;
}
