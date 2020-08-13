/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

void n_factor_insert(n_factor_t * factors, mp_limb_t p, ulong exp)
{
   ulong i;
   
   for (i = 0; i < factors->num; i++)
   {
      if (factors->p[i] == p) break;
   }
   
   if (i != factors->num) 
   {
      factors->exp[i] += exp;
   } else
   {
      factors->p[i] = p;
      factors->exp[i] = exp;
      factors->num++;
   }
}
