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
/****************************************************************************

   Copyright (C) 2009 William Hart

*****************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

int n_is_probabprime(mp_limb_t n)
{
   if (n <= 1UL) return 0;
   if (n == 2UL) return 1;
   if ((n & 1UL) == 0) return 0;

   if (n_is_perfect_power235(n)) return 0;
   
#if FLINT64
   if (n >= 10000000000000000UL) return n_is_probabprime_BPSW(n);
#endif

	mp_limb_t d = n - 1;
   unsigned int norm;

   count_trailing_zeros(norm, d);
   d >>= norm;

#if FLINT64
   if (n < 1122004669633UL)
#else
   if (n < 2147483648UL)
#endif  
   {
      if (n < FLINT_ODDPRIME_SMALL_CUTOFF)
         return n_is_oddprime_small(n);

      n_compute_primes(74000);
      if (n < flint_primes_cutoff)
         return n_is_oddprime_binary(n);
      
      double npre = n_precompute_inverse(n);
      
	   if (n < 9080191UL) 
      { 
         if (n_is_strong_probabprime_precomp(n, npre, 31UL, d) 
            && n_is_strong_probabprime_precomp(n, npre, 73UL, d)) return 1;
         else return 0;
      }

#if FLINT64
      if (n < 4759123141UL)
      {
#endif
         if (n_is_strong_probabprime_precomp(n, npre, 2UL, d) 
            && n_is_strong_probabprime_precomp(n, npre, 7UL, d) 
            && n_is_strong_probabprime_precomp(n, npre, 61UL, d)) return 1;
         else return 0;
#if FLINT64
      }
      
      if (n_is_strong_probabprime_precomp(n, npre, 2UL, d) 
         && n_is_strong_probabprime_precomp(n, npre, 13UL, d) 
         && n_is_strong_probabprime_precomp(n, npre, 23UL, d) 
         && n_is_strong_probabprime_precomp(n, npre, 1662803UL, d)) 
            if (n != 46856248255981UL) return 1;

      return 0;
#endif
   }

	mp_limb_t ninv = n_preinvert_limb(n);
   
   if (n_is_strong_probabprime2_preinv(n, ninv, 2UL, d) 
      && n_is_strong_probabprime2_preinv(n, ninv, 3UL, d) 
      && n_is_strong_probabprime2_preinv(n, ninv, 7UL, d) 
      && n_is_strong_probabprime2_preinv(n, ninv, 61UL, d) 
      && n_is_strong_probabprime2_preinv(n, ninv, 24251UL, d)) 
         if (n != 46856248255981UL) return 1;

   return 0;
}
