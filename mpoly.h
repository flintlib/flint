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
#include "ulong_extras.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef enum {
   ORD_LEX, ORD_DEGLEX, ORD_DEGREVLEX
} ordering_t;

#define MPOLY_NUM_ORDERINGS 3

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
      case ORD_DEGREVLEX:                                             \
         (deg) = (rev) = 1; break;                                    \
      default:                                                        \
         flint_throw(FLINT_ERROR, "Invalid ordering in fmpz_mpoly");  \
      }                                                               \
   } while (0)

#define masks_from_bits_ord(maskhi, masklo, bits, ord)                       \
   do {                                                                      \
      if ((ord) == ORD_DEGREVLEX)                                            \
      {                                                                      \
         (masklo) = -WORD(1);                                                \
         (maskhi) = (WORD(1) << ((bits)*((FLINT_BITS)/(bits) - 1))) - 1;     \
      } else                                                                 \
      {                                                                      \
         (masklo) = 0;                                                       \
         (maskhi) = 0;                                                       \
      }                                                                      \
   } while (0)


#define words_per_exp(nfields, bits)                                  \
   (((nfields) - 1)/(FLINT_BITS/(bits)) + 1)


/* Orderings *****************************************************************/

MPOLY_INLINE
ordering_t mpoly_ordering_randtest(flint_rand_t state)
{
   return (ordering_t) n_randint(state, MPOLY_NUM_ORDERINGS);
}

MPOLY_INLINE
int mpoly_ordering_isdeg(ordering_t ord)
{
   return ord == ORD_DEGLEX || ord == ORD_DEGREVLEX;
}

MPOLY_INLINE
int mpoly_ordering_isrev(ordering_t ord)
{
   return ord == ORD_DEGREVLEX;
}

MPOLY_INLINE
void mpoly_ordering_print(ordering_t ord)
{
   switch (ord)
   {
   case ORD_LEX:
      printf("lex");
      break;
   case ORD_DEGLEX:
      printf("deglex");
      break;
   case ORD_DEGREVLEX:
      printf("degrevlex");
      break;
   default:
      printf("Unknown ordering in mpoly_ordering_print.");
   }
}

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
int mpoly_monomial_overflows(ulong * exp2, slong N, ulong mask)
{
   slong i;
   for (i = 0; i < N; i++)
   {
      if ((exp2[i] & mask) != 0)
         return 1;
   }
   return 0;
}

MPOLY_INLINE
int mpoly_monomial_overflows1(ulong exp, ulong mask)
{
   return (exp & mask) != 0;
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
int mpoly_monomial_divides_test(const ulong * exp2,
                                       const ulong * exp3, slong N, ulong mask)
{
   slong i;
   for (i = 0; i < N; i++)
      if (((exp2[i] - exp3[i]) & mask) != 0)
         return 0;

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
int mpoly_monomial_lt(const ulong * exp2, const ulong * exp3,
                                           slong N, ulong maskhi, ulong masklo)
{
   slong i = 0;

      if (exp2[i] != exp3[i])
         return (exp3[i]^maskhi) < (exp2[i]^maskhi);

   for (i++; i < N; i++)
   {
      if (exp2[i] != exp3[i])
         return (exp3[i]^masklo) < (exp2[i]^masklo);
   }

   return 0;
}

MPOLY_INLINE
int mpoly_monomial_gt(const ulong * exp2, const ulong * exp3,
                                           slong N, ulong maskhi, ulong masklo)
{
   slong i = 0;

      if (exp2[i] != exp3[i])
         return (exp3[i]^maskhi) > (exp2[i]^maskhi);

   for (i++; i < N; i++)
   {
      if (exp2[i] != exp3[i])
         return (exp3[i]^masklo) > (exp2[i]^masklo);
   }

   return 0;
}

MPOLY_INLINE
int mpoly_monomial_cmp(const ulong * exp2, const ulong * exp3,
                                           slong N, ulong maskhi, ulong masklo)
{
   slong i = 0;

      if (exp2[i] != exp3[i])
      {
         if ((exp2[i]^maskhi) > (exp3[i]^maskhi))
            return 1;
         else
            return -1;
      }

   for (i++; i < N; i++)
   {
      if (exp2[i] != exp3[i])
      {
         if ((exp2[i]^masklo) > (exp3[i]^masklo))
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

MPOLY_INLINE
void mpoly_max_degrees_tight(slong * max_exp,
                             ulong * exps, slong len, slong * prods, slong num)
{
   slong i, j;
   
   for (j = 0; j < num; j++)
      max_exp[j] = 0;
	  
   for (i = 0; i < len; i++)
   {
      for (j = 0; j < num; j++)
	  {
	     slong d1 = (exps[i] % prods[j + 1])/prods[j];
      
	     if (d1 > max_exp[j])
		    max_exp[j] = d1;
	  }
   }
}

/* Monomial arrays ***********************************************************/

FLINT_DLL slong mpoly_exp_bits(const ulong * user_exp, slong nfields, int deg);

FLINT_DLL slong mpoly_optimize_bits(slong bits, slong nfields);

FLINT_DLL void   mpoly_pack_vec(ulong * exp1, const ulong * exp2, slong bits, slong nfields, slong len);

FLINT_DLL void mpoly_unpack_vec(ulong * exp1, const ulong * exp2, slong bits, slong nfields, slong len);



FLINT_DLL void mpoly_get_monomial(ulong * exps, const ulong * poly_exps,
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL void mpoly_set_monomial(ulong * exp1, const ulong * exp2,
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL void mpoly_unpack_monomials(ulong * exps1, slong bits1,
                       const ulong * exps2, slong bits2, slong len, slong num);

FLINT_DLL void mpoly_pack_monomials_tight(ulong * exp1,
                  const ulong * exp2, slong len, const slong * mults, 
                                           slong num, slong extra, slong bits);

FLINT_DLL void mpoly_unpack_monomials_tight(ulong * e1, ulong * e2, slong len,
                            slong * mults, slong num, slong extra, slong bits);

FLINT_DLL int mpoly_monomial_exists(slong * index, const ulong * poly_exps,
            const ulong * exp, slong len, slong N, ulong maskhi, ulong masklo);

FLINT_DLL void mpoly_max_degrees(ulong * max_degs, const ulong * poly_exps,
                                               slong len, slong bits, slong n);

FLINT_DLL void mpoly_search_monomials(
                slong ** e_ind, ulong * e, slong * e_score,
                slong * t1, slong * t2, slong *t3,
                slong lower, slong upper,
                const ulong * a, slong a_len, const ulong * b, slong b_len,
                                          slong N, ulong maskhi, ulong masklo);

MPOLY_INLINE
int mpoly_monomials_test(ulong * exps, slong len, slong N,
                                                    ulong maskhi, ulong masklo)
{
   slong i;

   for (i = 0; i + 1 < len; i++)
   {
      if (!mpoly_monomial_gt(exps + (i + 1)*N, exps + i*N, N, maskhi, masklo))
         return 0;
   }
   return 1;
}

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
void * _mpoly_heap_pop1(mpoly_heap1_s * heap, slong * heap_len, ulong maskhi)
{
   ulong exp;
   slong i, j, s = --(*heap_len);
   void * x = heap[1].next;

   i = 1;
   j = 2;

   while (j < s)
   {
      if ((heap[j].exp^maskhi) <= (heap[j + 1].exp^maskhi))
         j++;
      heap[i] = heap[j];
      i = j;
      j = HEAP_LEFT(j);
   }

   /* insert last element into heap[i] */
   exp = heap[s].exp;
   j = HEAP_PARENT(i);

   while (i > 1 && (exp^maskhi) > (heap[j].exp^maskhi))
   {
     heap[i] = heap[j];
     i = j;
     j = HEAP_PARENT(j);
   }

   heap[i] = heap[s];

   return x;
}

MPOLY_INLINE
void _mpoly_heap_insert1(mpoly_heap1_s * heap, ulong exp, void * x,
                              slong * next_loc, slong * heap_len, ulong maskhi)
{
   slong i = *heap_len, j, n = *heap_len;

   if (i != 1 && exp == heap[1].exp)
   {
      ((mpoly_heap_t *) x)->next = heap[1].next;
      heap[1].next = x;

      return;
   }

   if (*next_loc < *heap_len)
   {
      if (exp == heap[*next_loc].exp)
      {
         ((mpoly_heap_t *) x)->next = heap[*next_loc].next;
         heap[*next_loc].next = x;
         return;
      }
   }



   while ((j = HEAP_PARENT(i)) >= 1)
   {
      if (exp == heap[j].exp)
      {
         ((mpoly_heap_t *) x)->next = heap[j].next;
         heap[j].next = x;

         *next_loc = j;


         return;
      } else if ((exp^maskhi) > (heap[j].exp^maskhi))
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
void * _mpoly_heap_pop(mpoly_heap_s * heap, slong * heap_len, slong N,
                                                    ulong maskhi, ulong masklo)
{
   ulong * exp;
   slong i, j, s = --(*heap_len);
   mpoly_heap_t * x = heap[1].next;

   i = 1;
   j = 2;

   while (j < s)
   {
      if (!mpoly_monomial_gt(heap[j + 1].exp, heap[j].exp, N, maskhi, masklo))
         j++;
      heap[i] = heap[j];
      i = j;
      j = HEAP_LEFT(j);
   }

   /* insert last element into heap[i] */
   exp = heap[s].exp;
   j = HEAP_PARENT(i);

   while (i > 1 && mpoly_monomial_gt(heap[j].exp, exp, N, maskhi, masklo))
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
       slong * next_loc, slong * heap_len, slong N, ulong maskhi, ulong masklo)
{
   slong i = *heap_len, j, n = *heap_len;

   if (i != 1 && mpoly_monomial_equal(exp, heap[1].exp, N))
   {
      ((mpoly_heap_t *) x)->next = heap[1].next;
      heap[1].next = x;

      return 0;
   }

   if (*next_loc < *heap_len)
   {
      if (mpoly_monomial_equal(exp, heap[*next_loc].exp, N))
      {
         ((mpoly_heap_t *) x)->next = heap[*next_loc].next;
         heap[*next_loc].next = x;
         return 0;
      }
   }

   while ((j = HEAP_PARENT(i)) >= 1)
   {
      if (!mpoly_monomial_gt(heap[j].exp, exp, N, maskhi, masklo))
         break;

      i = j;
   }

   if (j >= 1 && mpoly_monomial_equal(exp, heap[j].exp, N))
   {
      ((mpoly_heap_t *) x)->next = heap[j].next;
      heap[j].next = x;
      *next_loc = j;

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

/* chunking */

/*
   Set i1[i] to the index of the i-th "coefficient" in variable k of num
   variables, each taking the given number of bits in the exponent. This
   assumes there are l1 "coefficients" in a list of len1 exponents.
   Note this doesn't currently mask the relevant bits.
*/ 
MPOLY_INLINE
void mpoly_main_variable_terms1(slong * i1, slong * n1, const ulong * exp1,
                          slong l1, slong len1, slong k, slong num, slong bits)
{
   slong i, j = 0;
   slong shift = bits*(FLINT_BITS/bits - (num - k + 1));

   i1[0] = 0;
   for (i = 0; i < l1 - 1; i++)
   {
      while (j < len1 && (l1 - i - 1) == (slong) (exp1[j] >> shift))
         j++;

      i1[i + 1] = j;
      n1[i] = j - i1[i];
   }
   n1[l1 - 1] = len1 - j;
}

#ifdef __cplusplus
}
#endif

#endif
