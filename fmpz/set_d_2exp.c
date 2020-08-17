/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <math.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_set_d_2exp(fmpz_t f, double m, slong exp)
{
   int exp2;
   
   m = frexp(m, &exp2);
   exp += exp2;

   if (exp >= 53)
   {
      fmpz_set_d(f, ldexp(m, 53));
      fmpz_mul_2exp(f, f, exp - 53);
   } else if (exp < 0)
      fmpz_set_ui(f, 0);
   else
      fmpz_set_d(f, ldexp(m, exp));    
}

