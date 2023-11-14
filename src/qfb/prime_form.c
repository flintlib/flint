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

void qfb_prime_form(qfb_t r, fmpz_t D, fmpz_t p)
{
   fmpz_t q, rem, s, t;

   fmpz_init(s);

   if (fmpz_cmp_ui(p, 2) == 0) /* special case, p = 2 */
   {
      ulong m8 = fmpz_fdiv_ui(D, 8);
      if (m8 == 4)
         fmpz_set_ui(r->b, 2);
      else
         fmpz_set_ui(r->b, m8);

      fmpz_sub_ui(s, D, m8);
      fmpz_neg(s, s);
      fmpz_fdiv_q_2exp(r->c, s, 3);
      fmpz_set(r->a, p);

      fmpz_clear(s);

      return;
   }

   fmpz_init(t);

   fmpz_mod(t, D, p);
   if (fmpz_is_zero(t)) /* special case, p | D */
   {
      fmpz_init(q);
      fmpz_init(rem);

      fmpz_fdiv_q(s, D, p); /* s = D/p */
      if (fmpz_is_zero(s))
         fmpz_set(t, s);
      else
         fmpz_sub(t, p, s); /* t = -s mod p */
      while (fmpz_fdiv_ui(t, 4) != 0) /* ensure t is divisible by 4 */
         fmpz_add(t, t, p);
      fmpz_add(t, t, s); /* t = first possible value for d/p + 4c */
      fmpz_fdiv_q(t, t, p);
      fmpz_sqrtrem(q, rem, t);
      if (!fmpz_is_zero(rem)) /* t/p is not a square */
      {
         if (fmpz_is_even(t)) /* b/p is even as t/p is */
            fmpz_add_ui(q, q, 1 + fmpz_is_even(q));
         else
            fmpz_add_ui(q, q, 1 + fmpz_is_odd(q));
      } /* q = b/p */

      fmpz_mul(r->b, q, p);
      fmpz_mul(q, q, q);
      fmpz_mul(q, q, p);
      fmpz_sub(q, q, s);
      fmpz_fdiv_q_2exp(r->c, q, 2);
      fmpz_set(r->a, p);

      fmpz_clear(q);
      fmpz_clear(rem);
   } else
   {
      fmpz_sqrtmod(t, t, p);
      fmpz_sub(s, D, t);

      if (fmpz_is_odd(s))
         fmpz_sub(t, p, t);

      fmpz_set(r->a, p);
      fmpz_set(r->b, t);
      fmpz_mul(r->c, r->b, r->b);
      fmpz_sub(r->c, r->c, D);
      fmpz_divexact(r->c, r->c, r->a);
      fmpz_fdiv_q_2exp(r->c, r->c, 2);
   }

   fmpz_clear(s);
   fmpz_clear(t);
}
