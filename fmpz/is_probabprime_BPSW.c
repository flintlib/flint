/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

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
