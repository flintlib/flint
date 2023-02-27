/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 William Hart
    Copyright (C) 2020 Chia Network Inc

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "qfb.h"

/*
   This is based (with permission) on an implementation of NUCOMP
   by Michael Jacobson which is an adaption to binary quadratic
   forms of the algorithm presented in Appendix A of "Solving the
   Pell Equation", Michael Jacobson and High Williams, CMS Books in
   Mathematics, Springer 2009.
*/
void qfb_nucomp(qfb_t r, const qfb_t f, const qfb_t g, fmpz_t D, fmpz_t L)
{
   fmpz_t a1, a2, c2, ca, cb, cc, k, s, sp, ss, m, t, u2, v1, v2;

   if (fmpz_cmp(f->a, g->a) > 0)
   {
      qfb_nucomp(r, g, f, D, L);
      return;
   }

   fmpz_init(a1); fmpz_init(a2); fmpz_init(c2);
   fmpz_init(ca); fmpz_init(cb); fmpz_init(cc);
   fmpz_init(k); fmpz_init(m);
   fmpz_init(s); fmpz_init(sp); fmpz_init(ss);
   fmpz_init(t); fmpz_init(u2); fmpz_init(v1); fmpz_init(v2);

   /* nucomp calculation */

   fmpz_set(a1, f->a);
   fmpz_set(a2, g->a);
   fmpz_set(c2, g->c);

   fmpz_add(ss, f->b, g->b);
   fmpz_fdiv_q_2exp(ss, ss, 1);

   fmpz_sub(m, f->b, g->b);
   fmpz_fdiv_q_2exp(m, m, 1);

   fmpz_fdiv_r(t, a2, a1);
   if (fmpz_is_zero(t))
   {
      fmpz_set_ui(v1, 0);
      fmpz_set(sp, a1);
   } else
      fmpz_gcdinv(sp, v1, t, a1);

   fmpz_mul(k, m, v1);
   fmpz_fdiv_r(k, k, a1);

   if (!fmpz_is_one(sp))
   {
      fmpz_xgcd(s, v2, u2, ss, sp);

      fmpz_mul(k, k, u2);
      fmpz_mul(t, v2, c2);
      fmpz_sub(k, k, t);

      if (!fmpz_is_one(s))
      {
         fmpz_divexact(a1, a1, s);
         fmpz_divexact(a2, a2, s);
         fmpz_mul(c2, c2, s);
      }

      fmpz_fdiv_r(k, k, a1);
   }

   if (fmpz_cmp(a1, L) < 0)
   {
      fmpz_mul(t, a2, k);

      fmpz_mul(ca, a2, a1);

      fmpz_mul_2exp(cb, t, 1);
      fmpz_add(cb, cb, g->b);

      fmpz_add(cc, g->b, t);
      fmpz_mul(cc, cc, k);
      fmpz_add(cc, cc, c2);

      fmpz_divexact(cc, cc, a1);
   } else
   {
      fmpz_t m1, m2, r1, r2, co1, co2, temp;

      fmpz_init(m1); fmpz_init(m2); fmpz_init(r1); fmpz_init(r2);
      fmpz_init(co1); fmpz_init(co2); fmpz_init(temp);

      fmpz_set(r2, a1);
      fmpz_set(r1, k);

      fmpz_xgcd_partial(co2, co1, r2, r1, L);

      fmpz_mul(t, a2, r1);
      fmpz_mul(m1, m, co1);
      fmpz_add(m1, m1, t);
      fmpz_divexact(m1, m1, a1);

      fmpz_mul(m2, ss, r1);
      fmpz_mul(temp, c2, co1);
      fmpz_sub(m2, m2, temp);
      fmpz_divexact(m2, m2, a1);

      fmpz_mul(ca, r1, m1);
      fmpz_mul(temp, co1, m2);
      if (fmpz_sgn(co1) < 0)
         fmpz_sub(ca, ca, temp);
      else
         fmpz_sub(ca, temp, ca);

      fmpz_mul(cb, ca, co2);
      fmpz_sub(cb, t, cb);
      fmpz_mul_2exp(cb, cb, 1);
      fmpz_divexact(cb, cb, co1);
      fmpz_sub(cb, cb, g->b);
      fmpz_mul_2exp(temp, ca, 1);
      fmpz_fdiv_r(cb, cb, temp);

      fmpz_mul(cc, cb, cb);
      fmpz_sub(cc, cc, D);
      fmpz_divexact(cc, cc, ca);
      fmpz_fdiv_q_2exp(cc, cc, 2);

      if (fmpz_sgn(ca) < 0)
      {
         fmpz_neg(ca, ca);
         fmpz_neg(cc, cc);
      }

      fmpz_clear(m1); fmpz_clear(m2); fmpz_clear(r1); fmpz_clear(r2);
      fmpz_clear(co1); fmpz_clear(co2); fmpz_clear(temp);
   }

   fmpz_set(r->a, ca);
   fmpz_set(r->b, cb);
   fmpz_set(r->c, cc);

   fmpz_clear(ca); fmpz_clear(cb); fmpz_clear(cc);
   fmpz_clear(k); fmpz_clear(m);
   fmpz_clear(s); fmpz_clear(sp); fmpz_clear(ss);
   fmpz_clear(t); fmpz_clear(u2); fmpz_clear(v1); fmpz_clear(v2);
   fmpz_clear(a1); fmpz_clear(a2); fmpz_clear(c2);
}
