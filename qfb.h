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

#ifndef QFB_H
#define QFB_H

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct qfb
{
    fmpz_t a;
    fmpz_t b;
    fmpz_t c;
} qfb;

typedef qfb qfb_t[1];

static __inline__
void qfb_init(qfb_t q)
{
   fmpz_init(q->a);
   fmpz_init(q->b);
   fmpz_init(q->c);
}

static __inline__
void qfb_clear(qfb_t q)
{
   fmpz_clear(q->a);
   fmpz_clear(q->b);
   fmpz_clear(q->c);
}

static __inline__
int qfb_equal(qfb_t f, qfb_t g)
{
   return (fmpz_equal(f->a, g->a) 
        && fmpz_equal(f->b, g->b) 
        && fmpz_equal(f->c, g->c));
}

static __inline__
void qfb_set(qfb_t f, qfb_t g)
{
   fmpz_set(f->a, g->a); 
   fmpz_set(f->b, g->b);
   fmpz_set(f->c, g->c);
}

static __inline__
void qfb_discriminant(fmpz_t D, qfb_t f)
{
   fmpz_t t;
   fmpz_init(t);

   fmpz_mul(t, f->a, f->c);
   fmpz_mul_2exp(t, t, 2);
   fmpz_mul(D, f->b, f->b);
   fmpz_sub(D, D, t);

   fmpz_clear(t);
}

static __inline__
void qfb_print(qfb_t q)
{
   printf("(");
   fmpz_print(q->a); printf(", ");
   fmpz_print(q->b); printf(", ");
   fmpz_print(q->c); printf(")");
}

static __inline__
void qfb_array_clear(qfb ** forms, long num)
{
   long k;

   for (k = 0; k < num; k++)
   {
      fmpz_clear((*forms)[k].a);
      fmpz_clear((*forms)[k].b);
      fmpz_clear((*forms)[k].c);
   }
   flint_free(*forms);
}

void qfb_reduce(qfb_t r, qfb_t f, fmpz_t D);

int qfb_is_reduced(qfb_t r);

long qfb_reduced_forms(qfb ** forms, long d);

long qfb_reduced_forms_large(qfb ** forms, long d);

void qfb_nucomp(qfb_t r, qfb_t f, qfb_t g, fmpz_t L);

void qfb_nudupl(qfb_t r, qfb_t f, fmpz_t L);

void qfb_pow_ui(qfb_t r, qfb_t f, fmpz_t D, ulong exp);

static __inline__
void qfb_inverse(qfb_t r, qfb_t f)
{
   qfb_set(r, f);
   
   if (fmpz_equal(f->a, f->c)
    || fmpz_equal(f->a, f->b))
    return;

   fmpz_neg(r->b, r->b);
}

static __inline__
int qfb_is_principal_form(qfb_t f, fmpz_t D)
{
   if (!fmpz_is_one(f->a))
      return 0;

   if (fmpz_is_odd(D)) /* D = 1 mod 4 */
      return fmpz_is_one(f->b);

   return fmpz_is_zero(f->b); /* D = 0 mod 4 */
}

static __inline__
void qfb_principal_form(qfb_t f, fmpz_t D)
{
   fmpz_set_ui(f->a, 1);
   
   if (fmpz_is_odd(D)) /* D = 1 mod 4 */
      fmpz_set_ui(f->b, 1);
   else /* D = 0 mod 4 */
      fmpz_set_ui(f->b, 0);

   fmpz_sub(f->c, f->b, D);
   fmpz_fdiv_q_2exp(f->c, f->c, 2);
}

static __inline__
int qfb_is_primitive(qfb_t f)
{
   fmpz_t g;
   int res;

   fmpz_init(g);

   fmpz_gcd(g, f->a, f->b);
   fmpz_gcd(g, g, f->c);
   res = fmpz_is_pm1(g);

   fmpz_clear(g);
   return res;
}

void qfb_prime_form(qfb_t r, fmpz_t D, fmpz_t p);

#ifdef __cplusplus
}
#endif

#endif
