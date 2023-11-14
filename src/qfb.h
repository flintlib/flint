/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef QFB_H
#define QFB_H

#ifdef QFB_INLINES_C
#define QFB_INLINE
#else
#define QFB_INLINE static inline
#endif

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

typedef struct
{
   qfb_t q;
   qfb_t q2;
   slong iter;
} qfb_hash_t;

QFB_INLINE
void qfb_init(qfb_t q)
{
   fmpz_init(q->a);
   fmpz_init(q->b);
   fmpz_init(q->c);
}

QFB_INLINE
void qfb_clear(qfb_t q)
{
   fmpz_clear(q->a);
   fmpz_clear(q->b);
   fmpz_clear(q->c);
}

QFB_INLINE
int qfb_equal(qfb_t f, qfb_t g)
{
   return (fmpz_equal(f->a, g->a)
        && fmpz_equal(f->b, g->b)
        && fmpz_equal(f->c, g->c));
}

QFB_INLINE
void qfb_set(qfb_t f, qfb_t g)
{
   fmpz_set(f->a, g->a);
   fmpz_set(f->b, g->b);
   fmpz_set(f->c, g->c);
}

QFB_INLINE
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

void qfb_print(qfb_t q);

QFB_INLINE
void qfb_array_clear(qfb ** forms, slong num)
{
   slong k;

   for (k = 0; k < num; k++)
   {
      fmpz_clear((*forms)[k].a);
      fmpz_clear((*forms)[k].b);
      fmpz_clear((*forms)[k].c);
   }
   flint_free(*forms);
}

qfb_hash_t * qfb_hash_init(slong depth);

void qfb_hash_clear(qfb_hash_t * qhash, slong depth);

void qfb_hash_insert(qfb_hash_t * qhash, qfb_t q,
                     qfb_t q2, slong iter, slong depth);

slong qfb_hash_find(qfb_hash_t * qhash, qfb_t q, slong depth);

void qfb_reduce(qfb_t r, qfb_t f, fmpz_t D);

int qfb_is_reduced(qfb_t r);

slong qfb_reduced_forms(qfb ** forms, slong d);

slong qfb_reduced_forms_large(qfb ** forms, slong d);

void qfb_nucomp(qfb_t r, const qfb_t f, const qfb_t g, fmpz_t D, fmpz_t L);

void qfb_nudupl(qfb_t r, const qfb_t f, fmpz_t D, fmpz_t L);

void qfb_pow_ui(qfb_t r, qfb_t f, fmpz_t D, ulong exp);

void qfb_pow(qfb_t r, qfb_t f, fmpz_t D, fmpz_t exp);

void qfb_pow_with_root(qfb_t r, qfb_t f, fmpz_t D, fmpz_t e, fmpz_t L);

QFB_INLINE
void qfb_inverse(qfb_t r, qfb_t f)
{
   qfb_set(r, f);

   if (fmpz_equal(f->a, f->c)
    || fmpz_equal(f->a, f->b))
    return;

   fmpz_neg(r->b, r->b);
}

QFB_INLINE
int qfb_is_principal_form(qfb_t f, fmpz_t D)
{
   if (!fmpz_is_one(f->a))
      return 0;

   if (fmpz_is_odd(D)) /* D = 1 mod 4 */
      return fmpz_is_one(f->b);

   return fmpz_is_zero(f->b); /* D = 0 mod 4 */
}

QFB_INLINE
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

QFB_INLINE
int qfb_is_primitive(qfb_t f)
{
   fmpz_t g;
   int res;

   fmpz_init(g);
   fmpz_gcd3(g, f->a, f->b, f->c);
   res = fmpz_is_pm1(g);
   fmpz_clear(g);

   return res;
}

void qfb_prime_form(qfb_t r, fmpz_t D, fmpz_t p);

int qfb_exponent_element(fmpz_t exponent, qfb_t f,
                                          fmpz_t n, ulong B1, ulong B2_sqrt);

int qfb_exponent(fmpz_t exponent, fmpz_t n, ulong B1, ulong B2_sqrt, slong c);

int qfb_exponent_grh(fmpz_t exponent, fmpz_t n, ulong B1, ulong B2_sqrt);

#ifdef __cplusplus
}
#endif

#endif
