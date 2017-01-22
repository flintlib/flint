/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qsieve.h"

void qsieve_poly_clear(qs_t qs_inf)
{
   slong i;

   fmpz_clear(qs_inf->A0);
   fmpz_clear(qs_inf->A);
   fmpz_clear(qs_inf->B);
   fmpz_clear(qs_inf->target_A);
   fmpz_clear(qs_inf->low_bound);
   fmpz_clear(qs_inf->upp_bound);

   for (i = 0; i < qs_inf->s; i++)
   {
       fmpz_clear(qs_inf->A0_divp[i]);
       fmpz_clear(qs_inf->B_terms[i]);
   }

   flint_free(qs_inf->B_terms);
   flint_free(qs_inf->A_ind);
   flint_free(qs_inf->A0_divp);
   flint_free(qs_inf->B0_terms);
   flint_free(qs_inf->A0_inv);
   flint_free(qs_inf->soln1);
   flint_free(qs_inf->soln2);
   flint_free(qs_inf->curr_subset);

   if (qs_inf->A_inv2B != NULL)
   {
       for (i = 0; i < qs_inf->s; i++)
       {
           flint_free(qs_inf->A_inv2B[i]);
       }
   }

   flint_free(qs_inf->A_inv2B);

#if HAVE_OPENMP
   for (i = 0; i < omp_get_max_threads(); i++)
   {
      fmpz_clear(qs_inf->poly[i].B);
      flint_free(qs_inf->poly[i].posn1);
      flint_free(qs_inf->poly[i].posn2);
      flint_free(qs_inf->poly[i].soln1);
      flint_free(qs_inf->poly[i].soln2);
      flint_free(qs_inf->poly[i].small);
      flint_free(qs_inf->poly[i].factor);
   }
#else
   fmpz_clear(qs_inf->poly[0].B);
   flint_free(qs_inf->poly[0].posn1);
   flint_free(qs_inf->poly[0].posn2);
   flint_free(qs_inf->poly[0].soln1);
   flint_free(qs_inf->poly[0].soln2);
   flint_free(qs_inf->poly[0].small);
   flint_free(qs_inf->poly[0].factor);
#endif

   qs_inf->B_terms = NULL;
   qs_inf->A_ind = NULL;
   qs_inf->A0_divp = NULL;
   qs_inf->B0_terms = NULL;
   qs_inf->A0_inv = NULL;
   qs_inf->soln1 = NULL;
   qs_inf->soln2 = NULL;
   qs_inf->A_inv2B = NULL;
   qs_inf->curr_subset = NULL;
}


