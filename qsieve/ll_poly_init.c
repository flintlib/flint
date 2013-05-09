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

    Copyright (C) 2006, 2011 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"

mp_limb_t qsieve_ll_poly_init(qs_t qs_inf)
{
   ulong num_primes = qs_inf->num_primes;
   ulong s = qs_inf->s; /* number of prime factors in A coeff */
   
   mp_limb_t ** A_inv2B;

   len_t i; 
        
   qs_inf->B_terms = flint_malloc(4*s*sizeof(mp_limb_t));
   qs_inf->A_ind = qs_inf->B_terms + s;  
   qs_inf->A_modp = qs_inf->A_ind + s;  
   qs_inf->inv_p2 = qs_inf->A_modp + s;

   qs_inf->A_inv2B = flint_malloc(s*sizeof(mp_limb_t *));

   qs_inf->A_inv = flint_malloc(3*num_primes*sizeof(mp_limb_t));  
   qs_inf->soln1 = qs_inf->A_inv + num_primes; 
   qs_inf->soln2 = qs_inf->soln1 + num_primes; 
   
   A_inv2B = qs_inf->A_inv2B;
   
   A_inv2B[0] = flint_malloc(num_primes*s*sizeof(mp_limb_t));
   for (i = 1; i < s; i++)
      A_inv2B[i] = A_inv2B[i - 1] + num_primes;
   
   return 0;
}
