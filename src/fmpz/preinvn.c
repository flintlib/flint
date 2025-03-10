/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"

void fmpz_preinvn_init(fmpz_preinvn_t inv, const fmpz_t f)
{
   fmpz c = *f;
   flint_bitcnt_t norm;
   nn_ptr t;

   if (c == 0)
   {
      flint_throw(FLINT_ERROR, "Exception (fmpz_preinvn_init). Division by zero.\n");
   } else if (!COEFF_IS_MPZ(c)) /* c is small */
   {
      ulong cc;
      inv->dinv = flint_malloc(sizeof(ulong));
      cc = FLINT_UABS(c);
      norm = flint_clz(cc);
      if (norm) cc <<= norm;
      flint_mpn_preinvn(inv->dinv, &cc, 1);
      inv->n = 1;
   } else /* c is big */
   {
      mpz_ptr mc = COEFF_TO_PTR(c);
      slong size = FLINT_ABS(mc->_mp_size);
      inv->dinv = flint_malloc(size*sizeof(ulong));
      norm = flint_clz(mc->_mp_d[size - 1]);
      if (norm)
      {
         t = flint_malloc(size*sizeof(ulong));
         mpn_lshift(t, mc->_mp_d, size, norm);
      } else
         t = mc->_mp_d;

      flint_mpn_preinvn(inv->dinv, t, size);

      inv->n = size;
      if (norm) flint_free(t);
   }

   inv->norm = norm;
}

void fmpz_preinvn_clear(fmpz_preinvn_t inv)
{
   flint_free(inv->dinv);
}
