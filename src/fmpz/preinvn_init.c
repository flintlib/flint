/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"

void fmpz_preinvn_init(fmpz_preinvn_t inv, const fmpz_t f)
{
   fmpz c = *f;
   flint_bitcnt_t norm;
   mp_ptr t;

   if (c == 0)
   {
      flint_throw(FLINT_ERROR, "Exception (fmpz_preinvn_init). Division by zero.\n");
   } else if (!COEFF_IS_MPZ(c)) /* c is small */
   {
      inv->dinv = flint_malloc(sizeof(mp_limb_t));
      if (c < 0) c = -c;
      norm = flint_clz(c);
      if (norm) c <<= norm;
      flint_mpn_preinvn(inv->dinv, (mp_ptr) &c, 1);
      inv->n = 1;
   } else /* c is big */
   {
      __mpz_struct * mc = COEFF_TO_PTR(c);
      slong size = FLINT_ABS(mc->_mp_size);
      inv->dinv = flint_malloc(size*sizeof(mp_limb_t));
      norm = flint_clz(mc->_mp_d[size - 1]);
      if (norm)
      {
         t = flint_malloc(size*sizeof(mp_limb_t));
         mpn_lshift(t, mc->_mp_d, size, norm);
      } else
         t = mc->_mp_d;

      flint_mpn_preinvn(inv->dinv, t, size);

      inv->n = size;
      if (norm) flint_free(t);
   }

   inv->norm = norm;
}
