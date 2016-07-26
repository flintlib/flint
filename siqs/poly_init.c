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
#include "siqs.h"

mp_limb_t qsieve_poly_init(qs_t qs_inf)
{
   ulong num_primes = qs_inf->num_primes;
   ulong s = qs_inf->s;    /* number of prime factors in A coeff */
   mp_limb_t ** A_inv2B;
   slong i;

   fmpz_init(qs_inf->A);
   fmpz_init(qs_inf->A0);
   fmpz_init(qs_inf->B);
   fmpz_init(qs_inf->C);
   fmpz_init(qs_inf->upp_bound);
   fmpz_init(qs_inf->low_bound);

   qs_inf->curr_subset = flint_malloc(s * sizeof(mp_limb_t));
   qs_inf->B_terms = flint_malloc(s * sizeof(mp_limb_t));
   qs_inf->A_ind = flint_malloc(s * sizeof(mp_limb_t));
   qs_inf->A0_divp = flint_malloc(s * sizeof(mp_limb_t));
   qs_inf->B0_terms = flint_malloc(s * sizeof(mp_limb_t));

   qs_inf->A_inv2B = flint_malloc(s * sizeof(mp_limb_t *));

   qs_inf->A0_inv = flint_malloc(num_primes * sizeof(mp_limb_t));
   qs_inf->soln1 = flint_malloc(num_primes * sizeof(mp_limb_t));
   qs_inf->soln2 = flint_malloc(num_primes * sizeof(mp_limb_t));
   qs_inf->posn1 = flint_malloc(num_primes * sizeof(mp_limb_t));
   qs_inf->posn2 = flint_malloc(num_primes * sizeof(mp_limb_t));

   A_inv2B = qs_inf->A_inv2B;

   for (i = 0; i < s; i++)
      A_inv2B[i] = flint_malloc(num_primes * sizeof(mp_limb_t));

   for (i = 0; i < s; i++)
   {
       fmpz_init(qs_inf->A0_divp[i]);
       fmpz_init(qs_inf->B_terms[i]);
   }

   return 0;
}

