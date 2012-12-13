/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 William Hart

******************************************************************************/

#undef ulong /* prevent clash with stdlib */
#include <stdlib.h>
#define ulong unsigned long
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qfb.h"

void qfb_pow(qfb_t r, qfb_t f, fmpz_t D, fmpz_t exp)
{
   fmpz_t L;
   qfb_t pow;

   if (fmpz_is_zero(exp))
   {
      qfb_principal_form(r, D);
      return;
   }

   if (fmpz_is_one(exp))
   {
      qfb_set(r, f);
      return;
   }

   fmpz_init(L);
   fmpz_abs(L, D);
   fmpz_root(L, L, 4);

   qfb_init(pow);
   
   qfb_set(pow, f);
   while (fmpz_is_even(exp))
   {
      qfb_nudupl(pow, pow, L);
      qfb_reduce(pow, pow, D);
      fmpz_fdiv_q_2exp(exp, exp, 1);
   }

   qfb_set(r, pow);
   fmpz_fdiv_q_2exp(exp, exp, 1);
   
   while (!fmpz_is_zero(exp))
   {
      qfb_nudupl(pow, pow, L);
      qfb_reduce(pow, pow, D);
      if (fmpz_is_odd(exp))
      {
         qfb_nucomp(r, r, pow, L);
         qfb_reduce(r, r, D);
      }
      fmpz_fdiv_q_2exp(exp, exp, 1);
   }

   qfb_clear(pow);
   fmpz_clear(L);
}
