/*
    Copyright (C) 2012 William Hart
    Copyright (C) 2020 Chia Network Inc

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

static void
qfb_nudupl_gcdinv(fmpz_t d, fmpz_t a, fmpz_t t, const fmpz_t f, const fmpz_t g)
{
   if (fmpz_cmp(f, g) < 0)
      fmpz_gcdinv(d, a, f, g);
   else
   {
      fmpz_fdiv_r(t, f, g);
      fmpz_gcdinv(d, a, t, g);
   }
}

void qfb_nudupl(qfb_t r, const qfb_t f, fmpz_t D, fmpz_t L)
{
   fmpz_t a1, b1, c1, ca, cb, cc, k, s, t, u2, v1, v2;

   fmpz_init(a1); fmpz_init(b1); fmpz_init(c1);
   fmpz_init(ca); fmpz_init(cb); fmpz_init(cc);
   fmpz_init(k);
   fmpz_init(s);
   fmpz_init(t); fmpz_init(u2); fmpz_init(v1); fmpz_init(v2);

   /* nucomp calculation */

   fmpz_set(a1, f->a);
   fmpz_set(c1, f->c);

   fmpz_zero(k);

   if (fmpz_cmpabs(b1, a1) == 0)
   {
      fmpz_set(s, a1);
      fmpz_zero(v2);
   } else if (fmpz_sgn(f->b) < 0)
   {
      fmpz_neg(b1, f->b);
      qfb_nudupl_gcdinv(s, v2, t, b1, a1);
      fmpz_neg(v2, v2);
   } else
      qfb_nudupl_gcdinv(s, v2, t, f->b, a1);

   fmpz_mul(t, v2, c1);
   fmpz_neg(k, t);

   if (!fmpz_is_one(s))
   {
      fmpz_divexact(a1, a1, s);
      fmpz_mul(c1, c1, s);
   }

   fmpz_fdiv_r(k, k, a1);

   if (fmpz_cmp(a1, L) < 0)
   {
      fmpz_mul(t, a1, k);

      fmpz_mul(ca, a1, a1);

      fmpz_mul_2exp(cb, t, 1);
      fmpz_add(cb, cb, f->b);

      fmpz_add(cc, f->b, t);
      fmpz_mul(cc, cc, k);
      fmpz_add(cc, cc, c1);

      fmpz_divexact(cc, cc, a1);
   } else
   {
      fmpz_t m2, r1, r2, co1, co2, temp;

      fmpz_init(m2); fmpz_init(r1); fmpz_init(r2);
      fmpz_init(co1); fmpz_init(co2); fmpz_init(temp);

      fmpz_set(r2, a1);
      fmpz_set(r1, k);

      fmpz_xgcd_partial(co2, co1, r2, r1, L);

      fmpz_mul(t, a1, r1);

      fmpz_mul(m2, f->b, r1);
      fmpz_mul(temp, c1, co1);
      fmpz_sub(m2, m2, temp);
      fmpz_divexact(m2, m2, a1);

      fmpz_mul(ca, r1, r1);
      fmpz_mul(temp, co1, m2);
      if (fmpz_sgn(co1) < 0)
         fmpz_sub(ca, ca, temp);
      else
         fmpz_sub(ca, temp, ca);

      fmpz_mul(cb, ca, co2);
      fmpz_sub(cb, t, cb);
      fmpz_mul_2exp(cb, cb, 1);
      fmpz_divexact(cb, cb, co1);
      fmpz_sub(cb, cb, f->b);
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

      fmpz_clear(m2); fmpz_clear(r1); fmpz_clear(r2);
      fmpz_clear(co1); fmpz_clear(co2); fmpz_clear(temp);
   }

   fmpz_set(r->a, ca);
   fmpz_set(r->b, cb);
   fmpz_set(r->c, cc);

   fmpz_clear(ca); fmpz_clear(cb); fmpz_clear(cc);
   fmpz_clear(k);
   fmpz_clear(s);
   fmpz_clear(t); fmpz_clear(u2); fmpz_clear(v2);
   fmpz_clear(a1); fmpz_clear(b1); fmpz_clear(c1);
}
