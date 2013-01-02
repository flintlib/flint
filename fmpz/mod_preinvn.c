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

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "mpn_extras.h"

void fmpz_mod2n_preinvn(fmpz_t r, const fmpz_t f, const fmpz_t m, mp_srcptr dinv)
{
   mp_ptr a, b, mp, fp;
   __mpz_struct * rem;
   mp_bitcnt_t norm;
   long mn = fmpz_size(m);
   long fn = fmpz_size(f);
   mp_limb_t ts[30], cy;
   
   rem = _fmpz_promote(r);
   mpz_realloc(rem, 2*mn);
   
   mp = COEFF_TO_PTR(*m)->_mp_d;
   fp = COEFF_TO_PTR(*f)->_mp_d;
   count_leading_zeros(norm, mp[mn - 1]);

   a = rem->_mp_d;
   mpn_zero(a + fn, 2*mn - fn);

   if (norm)
   {
      cy = mpn_lshift(a, fp, fn, norm);
      if (cy)
      {
         a[fn] = cy;
         fn++;
      }
      
      if (mn > 30)
         b = flint_malloc(mn*sizeof(mp_limb_t));
      else
         b = ts;
      mpn_lshift(b, mp, mn, norm);
   } else
   {
      flint_mpn_copyi(a, fp, fn);
      b = mp;
   }
   
   flint_mpn_rem2n_preinvn(a, b, mn, dinv);

   if (norm)
   {
      mpn_rshift(a, a, fn, norm);
      if (mn > 30)
         flint_free(b);
   }

   MPN_NORM(a, fn);
   rem->_mp_size = fn;
   _fmpz_demote_val(r);
}
