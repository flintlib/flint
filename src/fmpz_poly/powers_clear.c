/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_powers_clear(fmpz ** powers, slong len)
{
   slong i;
   for (i = 0; i < 2*len - 1; i++)
      _fmpz_vec_clear(powers[i], len - 1);

   flint_free(powers);
}

void
fmpz_poly_powers_clear(fmpz_poly_powers_precomp_t pinv)
{
   _fmpz_poly_powers_clear(pinv->powers, pinv->len);
}
