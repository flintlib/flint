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

void fmpz_preinvn(mp_ptr finv, const fmpz_t f)
{
   long fn = fmpz_size(f);
   mp_ptr fp, f1;
   mp_bitcnt_t norm;

   fp = COEFF_TO_PTR(*f)->_mp_d;
   count_leading_zeros(norm, fp[fn - 1]);

   if (norm)
   {
      f1 = flint_malloc(fn*sizeof(mp_limb_t));
      mpn_lshift(f1, fp, fn, norm);
   } else
      f1 = fp;
   
   flint_mpn_preinvn(finv, f1, fn);

   if (norm)
      flint_free(f1);
}
