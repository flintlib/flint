/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void 
_fmpq_poly_powers_clear(fmpq_poly_struct * powers, slong len)
{
   slong i;
   for (i = 0; i < 2*len - 1; i++)
      fmpq_poly_clear(powers + i);

   flint_free(powers);
}

void 
fmpq_poly_powers_clear(fmpq_poly_powers_precomp_t pinv)
{
   _fmpq_poly_powers_clear(pinv->powers, pinv->len);
}
