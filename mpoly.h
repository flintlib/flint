/*
    Copyright (C) 2016-2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef MPOLY_H
#define MPOLY_H

#ifdef MPOLY_INLINES_C
#define MPOLY_INLINE FLINT_DLL
#else
#define MPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif


typedef struct mpoly_heap_t
{
   ulong i;
   ulong j;
   struct mpoly_heap_t * next;
} mpoly_heap_t;

typedef struct mpoly_nheap_t
{
   ulong i;
   ulong j;
   struct mpoly_nheap_t * next;
   slong p;
} mpoly_nheap_t;

typedef struct mpoly_heap1_s
{
   ulong exp;
   void * next;
} mpoly_heap1_s;

typedef struct mpoly_heap_s
{
   ulong * exp;
   void * next;
} mpoly_heap_s;

/* Macros ********************************************************************/

#define degrev_from_ord(deg, rev, ord)                                \
   (deg) = (rev) = 0;                                                 \
   do {                                                               \
      switch (ord)                                                    \
      {                                                               \
      case ORD_LEX:                                                   \
         break;                                                       \
      case ORD_DEGLEX:                                                \
         (deg) = 1; break;                                            \
      case ORD_REVLEX:                                                \
         (rev) = 1; break;                                            \
      case ORD_DEGREVLEX:                                             \
         (deg) = (rev) = 1; break;                                    \
      default:                                                        \
         flint_throw(FLINT_ERROR, "Invalid ordering in fmpz_mpoly");  \
      }                                                               \
   } while (0)

/*  Monomials ****************************************************************/

MPOLY_INLINE
void mpoly_monomial_add(ulong * exp_ptr, const ulong * exp2,
                                                   const ulong * exp3, slong N)
{
   slong i;
   for (i = 0; i < N; i++)
      exp_ptr[i] = exp2[i] + exp3[i];
}

MPOLY_INLINE
void mpoly_monomial_sub_with_borrow(ulong * exp_ptr, const ulong * exp2, 
                                                   const ulong * exp3, slong N)
{
   slong i;
   ulong bw = 0;
   for (i = N - 1; i >= 0; i--)
      sub_ddmmss(bw, exp_ptr[i], WORD(0), exp2[i], WORD(0), exp3[i] - bw);  
}

MPOLY_INLINE
void mpoly_monomial_sub(ulong * exp_ptr, const ulong * exp2,
                                                   const ulong * exp3, slong N)
{
   slong i;
   for (i = 0; i < N; i++)
      exp_ptr[i] = exp2[i] - exp3[i];
}

MPOLY_INLINE
int mpoly_monomial_divides(ulong * exp_ptr, const ulong * exp2,
                                       const ulong * exp3, slong N, ulong mask)
{
   slong i;
   for (i = 0; i < N; i++)
   {
      exp_ptr[i] = exp2[i] - exp3[i];

      if ((exp_ptr[i] & mask) != 0)
         return 0;
   }

   return 1;
}

MPOLY_INLINE
int mpoly_monomial_divides1(ulong * exp_ptr, const ulong exp2,
                                                  const ulong exp3, ulong mask)
{
   (*exp_ptr) = exp2 - exp3;

   if (((exp2 - exp3) & mask) != 0)
      return 0;

   return 1;
}

MPOLY_INLINE
void mpoly_monomial_set(ulong * exp2, const ulong * exp3, slong N)
{
   slong i;
   for (i = 0; i < N; i++)
      exp2[i] = exp3[i];
}

MPOLY_INLINE
void mpoly_monomial_swap(ulong * exp2, ulong * exp3, slong N)
{
   slong i;
   ulong t;
   for (i = 0; i < N; i++)
   {
      t = exp2[i];
      exp2[i] = exp3[i];
      exp3[i] = t;
   }
}

MPOLY_INLINE
void mpoly_monomial_mul_si(ulong * exp2, const ulong * exp3, slong N, slong c)
{
   slong i;
   for (i = 0; i < N; i++)
      exp2[i] = exp3[i]*c;
}

MPOLY_INLINE
int mpoly_monomial_is_zero(const ulong * exp, slong N)
{
   slong i;
   for (i = 0; i < N; i++)
   {
      if (exp[i] != 0)
         return 0;
   }

   return 1;
}

MPOLY_INLINE
int mpoly_monomial_equal(const ulong * exp2, const ulong * exp3, slong N)
{
   slong i;

   for (i = 0; i < N; i++)
   {
      if (exp2[i] != exp3[i])
         return 0;
   }

   return 1;
}

MPOLY_INLINE
int mpoly_monomial_lt(const ulong * exp2, const ulong * exp3, slong N)
{
   slong i;

   for (i = 0; i < N; i++)
   {
      if (exp2[i] != exp3[i])
         return exp2[i] < exp3[i];
   }

   return 0;
}

MPOLY_INLINE
int mpoly_monomial_cmp(const ulong * exp2, const ulong * exp3, slong N)
{
   slong i;

   for (i = 0; i < N; i++)
   {
      if (exp2[i] != exp3[i])
      {
         if (exp2[i] > exp3[i])
            return 1;
         else
            return -1;
      }
   }

   return 0;
}

MPOLY_INLINE
int mpoly_monomial_divides_tight(slong e1, slong e2, slong * prods, slong num)
{
   slong j;

   for (j = 0; j < num; j++)
   {
      slong d1 = (e1 % prods[j + 1])/prods[j];
      slong d2 = (e2 % prods[j + 1])/prods[j];

      if (d1 < d2)
         return 0;
   }

   return 1;
}

/* Monomial arrays ***********************************************************/

FLINT_DLL void mpoly_get_monomial(ulong * exps, const ulong * poly_exps,
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL void mpoly_set_monomial(ulong * exp1, const ulong * exp2,
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL ulong * mpoly_unpack_monomials(slong bits1, const ulong * exps2, 
                                            slong len, slong num, slong bits2);

FLINT_DLL void mpoly_pack_monomials_tight(ulong * exp1,
                  const ulong * exp2, slong len, const slong * mults, 
                                           slong num, slong extra, slong bits);

FLINT_DLL void mpoly_unpack_monomials_tight(ulong * e1, ulong * e2, slong len,
                            slong * mults, slong num, slong extra, slong bits);

/* Heap **********************************************************************/

#define HEAP_LEFT(i) (2*(i))
#define HEAP_RIGHT(i) (2*(i) + 1)
#define HEAP_PARENT(i) ((i)/2)

#define HEAP_ASSIGN(h, c1, c2) \
   do {                                \
      (h).exp = (c1);                  \
      (h).next = (c2);                 \
   } while (0)

MPOLY_INLINE
void * _mpoly_heap_pop1(mpoly_heap1_s * heap, slong * heap_len)
{
   ulong exp;
   slong i, j, s = --(*heap_len);
   void * x = heap[1].next;

   i = 1;
   j = 2;

   while (j < s)
   {
      if (heap[j].exp >= heap[j + 1].exp)
         j++;
      heap[i] = heap[j];
      i = j;
      j = HEAP_LEFT(j);
   }

   /* insert last element into heap[i] */
   exp = heap[s].exp;
   j = HEAP_PARENT(i);

   while (i > 1 && exp < heap[j].exp)
   {
     heap[i] = heap[j];
     i = j;
     j = HEAP_PARENT(j);
   }

   heap[i] = heap[s];

   return x;
}

MPOLY_INLINE
void _mpoly_heap_insert1(mpoly_heap1_s * heap, ulong exp,
                                                    void * x, slong * heap_len)
{
   slong i = *heap_len, j, n = *heap_len;
   static slong next_loc = 0;

   if (i != 1 && exp == heap[1].exp)
   {
      ((mpoly_heap_t *) x)->next = heap[1].next;
      heap[1].next = x;

      return;
   }

   if (next_loc != 0 && next_loc < *heap_len)
   {
      if (exp == heap[next_loc].exp)
      {
         ((mpoly_heap_t *) x)->next = heap[next_loc].next;
         heap[next_loc].next = x;
         return;
      }
   }

   while ((j = HEAP_PARENT(i)) >= 1)
   {
      if (exp == heap[j].exp)
      {
         ((mpoly_heap_t *) x)->next = heap[j].next;
         heap[j].next = x;
         next_loc = j;

         return;
      } else if (exp < heap[j].exp)
         i = j;
      else
         break;
   }

   (*heap_len)++;

   while (n > i)
   {
      heap[n] = heap[HEAP_PARENT(n)];
      n = HEAP_PARENT(n);
   }

   HEAP_ASSIGN(heap[i], exp, x);
}

MPOLY_INLINE
void * _mpoly_heap_pop(mpoly_heap_s * heap, slong * heap_len, slong N)
{
   ulong * exp;
   slong i, j, s = --(*heap_len);
   mpoly_heap_t * x = heap[1].next;

   i = 1;
   j = 2;

   while (j < s)
   {
      if (mpoly_monomial_lt(heap[j + 1].exp, heap[j].exp, N))
         j++;
      heap[i] = heap[j];
      i = j;
      j = HEAP_LEFT(j);
   }

   /* insert last element into heap[i] */
   exp = heap[s].exp;
   j = HEAP_PARENT(i);

   while (i > 1 && mpoly_monomial_lt(exp, heap[j].exp, N))
   {
      heap[i] = heap[j];
      i = j;
      j = HEAP_PARENT(j);
   }

   heap[i] = heap[s];

   return x; 
}

MPOLY_INLINE
int _mpoly_heap_insert(mpoly_heap_s * heap, ulong * exp, void * x,
                                                     slong * heap_len, slong N)
{
   slong i = *heap_len, j, n = *heap_len;
   static slong next_loc = 0;

   if (i != 1 && mpoly_monomial_equal(exp, heap[1].exp, N))
   {
      ((mpoly_heap_t *) x)->next = heap[1].next;
      heap[1].next = x;

      return 0;
   }

   if (next_loc != 0 && next_loc < *heap_len)
   {
      if (mpoly_monomial_equal(exp, heap[next_loc].exp, N))
      {
         ((mpoly_heap_t *) x)->next = heap[next_loc].next;
         heap[next_loc].next = x;
         return 0;
      }
   }

   while ((j = HEAP_PARENT(i)) >= 1)
   {
      if (!mpoly_monomial_lt(exp, heap[j].exp, N))
         break;

      i = j;
   }

   if (j >= 1 && mpoly_monomial_equal(exp, heap[j].exp, N))
   {
      ((mpoly_heap_t *) x)->next = heap[j].next;
      heap[j].next = x;
      next_loc = j;

      return 0;
   }

   (*heap_len)++;

   while (n > i)
   {
      heap[n] = heap[HEAP_PARENT(n)];
      n = HEAP_PARENT(n);
   }

   HEAP_ASSIGN(heap[i], exp, x);

   return 1;
}

#ifdef __cplusplus
}
#endif

#endif
