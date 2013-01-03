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

void fmpz_mod_preinv1(fmpz_t r, const fmpz_t f, const fmpz_t m, mp_limb_t dinv)
{
   mp_limb_t * a, * b, * mp, * fp, cy;
   __mpz_struct * rem;
   mp_bitcnt_t norm;
   long mn = fmpz_size(m);
   long fn = fmpz_size(f);
   mp_limb_t t[71];
   
   if (fmpz_sgn(f) < 0 || fmpz_sgn(m) < 0
    || mn < 2 || fn > 70 || fn < mn)
   {
      fmpz_mod(r, f, m);
      return;
   }

   rem = _fmpz_promote(r);
   mpz_realloc(rem, fn + 1);
   
   mp = COEFF_TO_PTR(*m)->_mp_d;
   fp = COEFF_TO_PTR(*f)->_mp_d;
   count_leading_zeros(norm, mp[mn - 1]);

   a = rem->_mp_d;
   
   if (norm)
   {
      cy = mpn_lshift(a, fp, fn, norm);
      if (cy)
      {
         a[fn] = cy;
         fn++;
      }
      b = t + fn - mn;
      mpn_lshift(b, mp, mn, norm);
   } else
   {
      flint_mpn_copyi(a, fp, fn);
      b = mp;
   }
   
   flint_mpn_divrem_basecase_preinv1(t, a, fn, b, mn, dinv);

   if (norm)
      mpn_rshift(a, a, fn, norm);

   MPN_NORM(a, fn);
   rem->_mp_size = fn;
   _fmpz_demote_val(r);
}
