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

/*
   Shanks' NUCOMP as described in:
      "Computational aspects of NUCOMP", Michael J. Jacobson Jr.,
      Alfred J. van der Poorten, ANTS 2002, LNCS 2369, pp. 120--133.

   Computes the near reduced composition of forms f and g and returns
   the result in r.
*/
void qfb_nucomp(qfb_t r, qfb_t f, qfb_t g, fmpz_t L)
{
   int Fdivs;
   
   fmpz_t s, m, a, b, c, F, G, H, 
      Ax, Bx, By, Cy, Dy, x, y, 
      z, l, t1, t2, bx, by, 
      ax, ay, q, t, Q1, Q2, Q3, Q4,
      dx, dy, cx, cy;
   
   fmpz_init(s); fmpz_init(l); fmpz_init(m);
   fmpz_init(a); fmpz_init(c); fmpz_init(b);
   fmpz_init(F); fmpz_init(G); fmpz_init(H);
   fmpz_init(Ax); fmpz_init(Bx); fmpz_init(By); 
   fmpz_init(Cy); fmpz_init(Dy);
   fmpz_init(x); fmpz_init(y); fmpz_init(z);  
   fmpz_init(t1); fmpz_init(t2);
   fmpz_init(ax); fmpz_init(ay); fmpz_init(bx); fmpz_init(by); 
   fmpz_init(cx); fmpz_init(cy); fmpz_init(dx); fmpz_init(dy);
   fmpz_init(q); fmpz_init(t);
   fmpz_init(Q1); fmpz_init(Q2); fmpz_init(Q3); fmpz_init(Q4);
   
   /* Step 1: */
   if (fmpz_cmp(f->c, g->c) < 0)
   {
      qfb * temp = f;
      f = g;
      g = temp;
   }

   fmpz_add(s, f->b, g->b);
   fmpz_fdiv_q_2exp(s, s, 1);
   fmpz_sub(m, g->b, s);

   /* Step 2: */
   fmpz_xgcd(F, b, c, g->a, f->a);
   if ((Fdivs = fmpz_divisible(s, F)))
   {
      fmpz_set(G, F);
      fmpz_set(Ax, G);
      fmpz_mul(Bx, m, b);
   } else
   {
      /* Step 3: */
      fmpz_mod(t1, s, F);
      fmpz_gcdinv(G, y, t1, F);
      fmpz_divexact(H, F, G);
   }

   fmpz_divexact(By, f->a, G);
   fmpz_divexact(Cy, g->a, G);
   fmpz_divexact(Dy, s, G);
   
   /* Step 4 */
   if (!Fdivs)
   {
      fmpz_mod(t1, f->c, H);
      fmpz_mul(t1, t1, b);
      fmpz_mod(t2, g->c, H);
      fmpz_mul(t2, t2, c);
      fmpz_add(t1, t1, t2);
      fmpz_mod(t1, t1, H);
      fmpz_mul(t1, t1, y);
      fmpz_mod(l, t1, H);

      fmpz_divexact(t1, m, H);
      fmpz_mul(t1, t1, b);
      fmpz_divexact(t2, By, H);
      fmpz_mul(t2, t2, l);
      fmpz_add(Bx, t1, t2);
   }

   /* Step 5: */
   fmpz_mod(bx, Bx, By);
   fmpz_set(by, By);
   
   fmpz_set_ui(x, 1);
   fmpz_set_ui(y, 0);
   fmpz_set_ui(z, 0);

   while (1)
   {
      fmpz_abs(t, by);
      if (fmpz_cmp(t, L) <= 0 || fmpz_is_zero(bx))
      {
         if (fmpz_is_odd(z))
         {
            fmpz_neg(by, by);
            fmpz_neg(y, y);
         }

         fmpz_mul(ax, G, x);
         fmpz_mul(ay, G, y);

         break;
      }

      fmpz_fdiv_qr(q, t, by, bx);
      fmpz_set(by, bx);
      fmpz_set(bx, t);
      fmpz_mul(t, q, x);
      fmpz_sub(t, y, t);
      fmpz_set(y, x);
      fmpz_set(x, t);
      fmpz_add_ui(z, z, 1);
   }
   
   if (fmpz_is_zero(z))
   {
      /* Step 6: */
      fmpz_mul(Q1, Cy, bx);
      fmpz_sub(cx, Q1, m);
      fmpz_divexact(cx, cx, By);
      fmpz_mul(dx, bx, Dy);
      fmpz_sub(dx, dx, g->c);
      fmpz_divexact(dx, dx, By);
      fmpz_mul(r->a, by, Cy);
      fmpz_mul(t1, bx, cx);
      fmpz_mul(t2, G, dx);
      fmpz_sub(r->c, t1, t2);
      fmpz_set(t1, g->b);
      fmpz_mul_2exp(r->b, Q1, 1);
      fmpz_sub(r->b, t1, r->b);
   } else
   {
      /* Step 7: */
      fmpz_mul(t1, Cy, bx);
      fmpz_mul(t2, m, x);
      fmpz_sub(t1, t1, t2);
      fmpz_divexact(cx, t1, By);
      fmpz_mul(Q1, by, cx);
      fmpz_add(Q2, Q1, m);
      fmpz_mul(t1, Dy, bx);
      fmpz_mul(t2, g->c, x);
      fmpz_sub(t1, t1, t2);
      fmpz_divexact(dx, t1, By);
      fmpz_mul(Q3, y, dx);
      fmpz_add(Q4, Q3, Dy);
      fmpz_divexact(dy, Q4, x);
      if (!fmpz_is_zero(bx))
         fmpz_divexact(cy, Q2, bx);
      else
      {
         fmpz_mul(cy, cx, dy);
         fmpz_sub(cy, cy, f->c);
         fmpz_divexact(cy, cy, dx);
      }
      fmpz_mul(t1, by, cy);
      fmpz_mul(t2, ay, dy);
      fmpz_sub(r->a, t1, t2);
      fmpz_mul(t1, bx, cx);
      fmpz_mul(t2, ax, dx);
      fmpz_sub(r->c, t1, t2);
      fmpz_add(r->b, Q3, Q4);
      fmpz_mul(r->b, r->b, G);
      fmpz_sub(r->b, r->b, Q1);
      fmpz_sub(r->b, r->b, Q2);
   }

   fmpz_clear(s); fmpz_clear(l); fmpz_clear(m);
   fmpz_clear(a); fmpz_clear(b); fmpz_clear(c);
   fmpz_clear(F); fmpz_clear(G); fmpz_clear(H);
   fmpz_clear(Ax); fmpz_clear(Bx); fmpz_clear(By); 
   fmpz_clear(Cy); fmpz_clear(Dy);
   fmpz_clear(x); fmpz_clear(y); fmpz_clear(z);  
   fmpz_clear(t1); fmpz_clear(t2);
   fmpz_clear(ax); fmpz_clear(ay); fmpz_clear(bx); fmpz_clear(by); 
   fmpz_clear(cx); fmpz_clear(cy); fmpz_clear(dx); fmpz_clear(dy);
   fmpz_clear(q); fmpz_clear(t);
   fmpz_clear(Q1); fmpz_clear(Q2); fmpz_clear(Q3); fmpz_clear(Q4);
}
