/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"

int n_divides(ulong * q, ulong n, ulong p)
{
   ulong quo, rem;

   if (p == 0)
   {
      *q = 0;
      return n == 0;
   }

   /* purportedly the compiler optimises this */
   quo = n/p;
   rem = n%p;
   if (rem == 0)
   {
      *q = quo;
      return 1;
   } else
   {
      *q = 0;
      return 0;
   }
}
