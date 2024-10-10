/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "qsieve.h"

ulong qsieve_poly_init(qs_t qs_inf)
{
   ulong num_primes = qs_inf->num_primes;
   ulong s = qs_inf->s;    /* number of prime factors in A coeff */
   ulong ** A_inv2B;
   slong i;

   fmpz_init(qs_inf->A);
   fmpz_init(qs_inf->B);
   fmpz_init(qs_inf->upp_bound);
   fmpz_init(qs_inf->low_bound);

   qs_inf->curr_subset = flint_malloc(s * sizeof(ulong));
   qs_inf->first_subset = flint_malloc(s * sizeof(ulong));
   qs_inf->B_terms = flint_malloc(s * sizeof(ulong));
   qs_inf->A_ind = flint_malloc(s * sizeof(ulong));
   qs_inf->A_divp = flint_malloc(s * sizeof(ulong));
   qs_inf->B0_terms = flint_malloc(s * sizeof(ulong));

   qs_inf->A_inv2B = flint_malloc(s * sizeof(ulong *));

   qs_inf->A_inv = flint_malloc(num_primes * sizeof(ulong));
   qs_inf->soln1 = flint_malloc(num_primes * sizeof(ulong));
   qs_inf->soln2 = flint_malloc(num_primes * sizeof(ulong));

   qs_inf->poly = flint_malloc((qs_inf->num_handles + 1)* sizeof(qs_poly_s));

   for (i = 0; i <= qs_inf->num_handles ; i++)
   {
      fmpz_init(qs_inf->poly[i].B);
      qs_inf->poly[i].posn1 = flint_malloc((num_primes + 16)*sizeof(ulong));
      qs_inf->poly[i].posn2 = flint_malloc((num_primes + 16)*sizeof(ulong));
      qs_inf->poly[i].soln1 = flint_malloc((num_primes + 16)*sizeof(ulong));
      qs_inf->poly[i].soln2 = flint_malloc((num_primes + 16)*sizeof(ulong));
      qs_inf->poly[i].small = flint_malloc(qs_inf->small_primes*sizeof(ulong));
      qs_inf->poly[i].factor = flint_malloc(qs_inf->max_factors*sizeof(fac_t));
   }

   A_inv2B = qs_inf->A_inv2B;

   for (i = 0; (ulong) i < s; i++)
      A_inv2B[i] = flint_malloc(num_primes * sizeof(ulong));

   for (i = 0; (ulong) i < s; i++)
   {
       fmpz_init(qs_inf->A_divp[i]);
       fmpz_init(qs_inf->B_terms[i]);
   }

   return 0;
}

void qsieve_poly_clear(qs_t qs_inf)
{
   slong i;

   fmpz_clear(qs_inf->A);
   fmpz_clear(qs_inf->B);
   fmpz_clear(qs_inf->target_A);
   fmpz_clear(qs_inf->low_bound);
   fmpz_clear(qs_inf->upp_bound);

   for (i = 0; i < qs_inf->s; i++)
   {
       fmpz_clear(qs_inf->A_divp[i]);
       fmpz_clear(qs_inf->B_terms[i]);
   }

   flint_free(qs_inf->B_terms);
   flint_free(qs_inf->A_ind);
   flint_free(qs_inf->A_divp);
   flint_free(qs_inf->B0_terms);
   flint_free(qs_inf->A_inv);
   flint_free(qs_inf->soln1);
   flint_free(qs_inf->soln2);
   flint_free(qs_inf->curr_subset);
   flint_free(qs_inf->first_subset);

   if (qs_inf->A_inv2B != NULL)
   {
       for (i = 0; i < qs_inf->s; i++)
       {
           flint_free(qs_inf->A_inv2B[i]);
       }
   }

   flint_free(qs_inf->A_inv2B);

   for (i = 0; i <= qs_inf->num_handles; i++)
   {
      fmpz_clear(qs_inf->poly[i].B);
      flint_free(qs_inf->poly[i].posn1);
      flint_free(qs_inf->poly[i].posn2);
      flint_free(qs_inf->poly[i].soln1);
      flint_free(qs_inf->poly[i].soln2);
      flint_free(qs_inf->poly[i].small);
      flint_free(qs_inf->poly[i].factor);
   }
   flint_free(qs_inf->poly);

   qs_inf->B_terms = NULL;
   qs_inf->A_ind = NULL;
   qs_inf->A_divp = NULL;
   qs_inf->B0_terms = NULL;
   qs_inf->A_inv = NULL;
   qs_inf->soln1 = NULL;
   qs_inf->soln2 = NULL;
   qs_inf->A_inv2B = NULL;
   qs_inf->curr_subset = NULL;
}
