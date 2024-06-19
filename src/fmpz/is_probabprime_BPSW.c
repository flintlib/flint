/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

int fmpz_is_probabprime_BPSW(const fmpz_t n)
{
   fmpz_t b;
   int res = 1;

   fmpz_init_set_ui(b, 2);

   if (!fmpz_is_strong_probabprime(n, b) || !fmpz_is_probabprime_lucas(n))
      res = 0;

   fmpz_clear(b);

   return res;
}
