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

void qfb_pow(qfb_t r, qfb_t f, fmpz_t D, fmpz_t e)
{
   fmpz_t L, exp;
   qfb_t pow;

   if (fmpz_is_zero(e))
   {
      qfb_principal_form(r, D);
      return;
   }

   if (fmpz_is_one(e))
   {
      qfb_set(r, f);
      return;
   }

   fmpz_init(L);
   fmpz_init(exp);

   fmpz_set(exp, e);

   fmpz_abs(L, D);
   fmpz_root(L, L, 4);

   qfb_init(pow);

   qfb_set(pow, f);
   while (fmpz_is_even(exp))
   {
      qfb_nudupl(pow, pow, D, L);
      qfb_reduce(pow, pow, D);
      fmpz_fdiv_q_2exp(exp, exp, 1);
   }

   qfb_set(r, pow);
   fmpz_fdiv_q_2exp(exp, exp, 1);

   while (!fmpz_is_zero(exp))
   {
      qfb_nudupl(pow, pow, D, L);
      qfb_reduce(pow, pow, D);
      if (fmpz_is_odd(exp))
      {
         qfb_nucomp(r, r, pow, D, L);
         qfb_reduce(r, r, D);
      }
      fmpz_fdiv_q_2exp(exp, exp, 1);
   }

   qfb_clear(pow);
   fmpz_clear(L);
   fmpz_clear(exp);
}

void qfb_pow_with_root(qfb_t r, qfb_t f, fmpz_t D, fmpz_t e, fmpz_t L)
{
   fmpz_t exp;
   qfb_t pow;

   if (fmpz_is_zero(e))
   {
      qfb_principal_form(r, D);
      return;
   }

   if (fmpz_is_one(e))
   {
      qfb_set(r, f);
      return;
   }

   fmpz_init(exp);

   fmpz_set(exp, e);

   qfb_init(pow);

   qfb_set(pow, f);
   while (fmpz_is_even(exp))
   {
      qfb_nudupl(pow, pow, D, L);
      qfb_reduce(pow, pow, D);
      fmpz_fdiv_q_2exp(exp, exp, 1);
   }

   qfb_set(r, pow);
   fmpz_fdiv_q_2exp(exp, exp, 1);

   while (!fmpz_is_zero(exp))
   {
      qfb_nudupl(pow, pow, D, L);
      qfb_reduce(pow, pow, D);
      if (fmpz_is_odd(exp))
      {
         qfb_nucomp(r, r, pow, D, L);
         qfb_reduce(r, r, D);
      }
      fmpz_fdiv_q_2exp(exp, exp, 1);
   }

   qfb_clear(pow);
   fmpz_clear(exp);
}
