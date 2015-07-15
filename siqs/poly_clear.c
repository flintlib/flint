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

    Copyright (C) 2015 Nitin Kumar

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"

void qsieve_poly_clear(qs_t qs_inf)
{
   slong i;

   fmpz_clear(qs_inf->A0);
   fmpz_clear(qs_inf->A);
   fmpz_clear(qs_inf->B);
   fmpz_clear(qs_inf->C);
   fmpz_clear(qs_inf->target_A);

   for (i = 0; i < qs_inf->s; i++)
   {
       fmpz_clear(qs_inf->A0_divp[i]);
       fmpz_clear(qs_inf->B_terms[i]);
   }

   flint_free(qs_inf->q0_values);
   flint_free(qs_inf->B_terms);
   flint_free(qs_inf->A_ind);
   flint_free(qs_inf->A0_divp);
   flint_free(qs_inf->B0_terms);
   flint_free(qs_inf->A0_inv);
   flint_free(qs_inf->soln1);
   flint_free(qs_inf->soln2);
   flint_free(qs_inf->current_subset);

   if (qs_inf->A_inv2B != NULL)
   {
       for (i = 0; i < qs_inf->s; i++)
       {
           flint_free(qs_inf->A_inv2B[i]);
       }
   }

   flint_free(qs_inf->A_inv2B);

   qs_inf->q0_values = NULL;
   qs_inf->B_terms = NULL;
   qs_inf->A_ind = NULL;
   qs_inf->A0_divp = NULL;
   qs_inf->B0_terms = NULL;
   qs_inf->A0_inv = NULL;
   qs_inf->soln1 = NULL;
   qs_inf->soln2 = NULL;
   qs_inf->A_inv2B = NULL;
   qs_inf->current_subset = NULL;

}


