/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "qfb.h"

void qfb_pow_ui(qfb_t r, qfb_t f, fmpz_t D, ulong exp)
{
   fmpz_t L;
   qfb_t pow;

   if (exp == 0)
   {
      qfb_principal_form(r, D);
      return;
   }

   if (exp == 1)
   {
      qfb_set(r, f);
      return;
   }

   fmpz_init(L);
   fmpz_abs(L, D);
   fmpz_root(L, L, 4);

   qfb_init(pow);

   qfb_set(pow, f);
   while ((exp & 1) == 0)
   {
      qfb_nudupl(pow, pow, D, L);
      qfb_reduce(pow, pow, D);
      exp >>= 1;
   }

   qfb_set(r, pow);
   exp >>= 1;

   while (exp)
   {
      qfb_nudupl(pow, pow, D, L);
      qfb_reduce(pow, pow, D);
      if (exp & 1)
      {
         qfb_nucomp(r, r, pow, D, L);
         qfb_reduce(r, r, D);
      }
      exp >>= 1;
   }

   qfb_clear(pow);
   fmpz_clear(L);
}
