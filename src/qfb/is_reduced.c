/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qfb.h"

int qfb_is_reduced(qfb_t r)
{
   if (fmpz_cmp(r->c, r->a) < 0)
      return 0;

   if (fmpz_cmpabs(r->b, r->a) > 0)
      return 0;

   if (fmpz_cmpabs(r->a, r->b) == 0 || fmpz_cmp(r->a, r->c) == 0)
      if (fmpz_sgn(r->b) < 0)
         return 0;

   return 1;
}
