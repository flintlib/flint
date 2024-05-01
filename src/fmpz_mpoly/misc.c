/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mpoly.h"

void _fmpz_mpoly_sub_uiuiui_fmpz(ulong * c, const fmpz_t d)
{
   fmpz fc = *d;

   if (!COEFF_IS_MPZ(fc))
   {
        ulong f0, f1, f2;
        f0 = fc;
        f1 = f2 = FLINT_SIGN_EXT(f0);
        sub_dddmmmsss(c[2], c[1], c[0], c[2], c[1], c[0], f2, f1, f0);
   } else
   {
      slong size = fmpz_size(d);
      mpz_ptr m = COEFF_TO_PTR(fc);
      if (fmpz_sgn(d) < 0)
         mpn_add(c, c, 3, m->_mp_d, size);
      else
         mpn_sub(c, c, 3, m->_mp_d, size);
   }
}

void _fmpz_mpoly_add_uiuiui_fmpz(ulong * c, const fmpz_t d)
{
    fmpz fc = *d;

    if (!COEFF_IS_MPZ(fc))
    {
        ulong f0, f1, f2;
        f0 = fc;
        f1 = f2 = FLINT_SIGN_EXT(f0);
        add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], f2, f1, f0);
   } else
   {
      slong size = fmpz_size(d);
      mpz_ptr m = COEFF_TO_PTR(fc);
      if (fmpz_sgn(d) < 0)
         mpn_sub(c, c, 3, m->_mp_d, size);
      else
         mpn_add(c, c, 3, m->_mp_d, size);
   }
}
