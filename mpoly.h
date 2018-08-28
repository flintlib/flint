/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017 Daniel Schultz

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
#include "fmpz.h"
#include "ulong_extras.h"

#ifdef __cplusplus
 extern "C" {
#endif


#define MPOLY_MIN_BITS (UWORD(8))    /* minimum number of bits to pack into */

typedef enum {
   ORD_LEX, ORD_DEGLEX, ORD_DEGREVLEX
} ordering_t;

#define MPOLY_NUM_ORDERINGS 3

/* context *******************************************************************/

typedef struct
{
    slong nvars;    /* number of variables */
    slong nfields;  /* number of fields in exponent vector */
    ordering_t ord; /* monomial ordering */
    int deg;        /* is ord a degree ordering? */
    int rev;        /* is ord a reversed ordering? */
} mpoly_ctx_struct;

typedef mpoly_ctx_struct mpoly_ctx_t[1];

FLINT_DLL void mpoly_ctx_init(mpoly_ctx_t ctx, slong nvars, const ordering_t ord);

FLINT_DLL void mpoly_ctx_init_rand(mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars);

FLINT_DLL void mpoly_monomial_randbits_fmpz(fmpz * exp, flint_rand_t state, mp_bitcnt_t exp_bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_ctx_clear(mpoly_ctx_t mctx);

MPOLY_INLINE slong mpoly_words_per_exp(slong bits, const mpoly_ctx_t mctx)
{
    if (bits <= FLINT_BITS)
        return ((mctx->nfields) - 1)/(FLINT_BITS/(bits)) + 1;
    else
        return (bits + FLINT_BITS - 1)/FLINT_BITS*mctx->nfields;
}

/* heaps *********************************************************************/
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

/* trees *********************************************************************/
typedef struct mpoly_rbnode
{
    struct mpoly_rbnode * up;
    struct mpoly_rbnode * left;
    struct mpoly_rbnode * right;
    void * data;
    void * data2;
    slong key;
    int col;
} mpoly_rbnode_struct;

typedef mpoly_rbnode_struct mpoly_rbnode_t[1];

typedef struct mpoly_rbtree
{
    slong size;
    mpoly_rbnode_t head;  /* dummy node for pointer to head */
    mpoly_rbnode_t null;  /* dummy node to be pointed to by leaves */
} mpoly_rbtree_struct;

typedef mpoly_rbtree_struct mpoly_rbtree_t[1];

FLINT_DLL void mpoly_rbtree_init(mpoly_rbtree_t tree);

FLINT_DLL void mpoly_rbnode_clear(mpoly_rbtree_t tree, mpoly_rbnode_t node,
                                void ** dataout, slong * keysout, slong * idx);

FLINT_DLL void mpoly_rbtree_clear(mpoly_rbtree_t tree, void ** dataout, slong * keysout);

FLINT_DLL mpoly_rbnode_struct * mpoly_rbtree_get(int * new,
                                         struct mpoly_rbtree *tree, slong rcx);


/* Orderings *****************************************************************/

MPOLY_INLINE
ordering_t mpoly_ordering_randtest(flint_rand_t state)
{
   return (ordering_t) n_randint(state, MPOLY_NUM_ORDERINGS);
}

MPOLY_INLINE
int mpoly_ordering_isdeg(const mpoly_ctx_t mctx)
{
   return mctx->ord == ORD_DEGLEX || mctx->ord == ORD_DEGREVLEX;
}

MPOLY_INLINE
int mpoly_ordering_isrev(const mpoly_ctx_t mctx)
{
   return mctx->ord == ORD_DEGREVLEX;
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

/* Misc **********************************************************************/

FLINT_DLL void mpoly_gen_offset_shift(slong * _offset, slong * _shift,
                       slong idx, slong N, slong bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_gen_oneexp_offset_shift(ulong * oneexp, slong * offset, slong * shift,
                       slong var, slong N, slong bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_gen_oneexp_offset_mp(ulong * oneexp, slong * offset,
                       slong idx, slong N, slong bits, const mpoly_ctx_t mctx);

/*  Monomials ****************************************************************/

MPOLY_INLINE
void mpoly_monomial_zero(ulong * exp_ptr, slong N)
{
   slong i;
   for (i = 0; i < N; i++)
      exp_ptr[i] = 0;
}

MPOLY_INLINE
void mpoly_monomial_add(ulong * exp_ptr, const ulong * exp2,
                                                   const ulong * exp3, slong N)
{
   slong i;
   for (i = 0; i < N; i++)
      exp_ptr[i] = exp2[i] + exp3[i];
}

MPOLY_INLINE
void mpoly_monomial_add_mp(ulong * exp_ptr, const ulong * exp2,
                                                   const ulong * exp3, slong N)
{
    mpn_add_n(exp_ptr, exp2, exp3, N);
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
void mpoly_monomial_sub_mp(ulong * exp_ptr, const ulong * exp2,
                                                   const ulong * exp3, slong N)
{
    mpn_sub_n(exp_ptr, exp2, exp3, N);
}

MPOLY_INLINE
void mpoly_monomial_madd(ulong * exp1, const ulong * exp2, ulong scalar,
                                                   const ulong * exp3, slong N)
{
   slong i;
   for (i = 0; i < N; i++)
      exp1[i] = exp2[i] + scalar*exp3[i];
}

MPOLY_INLINE
void mpoly_monomial_madd_mp(ulong * exp1, const ulong * exp2, ulong scalar,
                                                   const ulong * exp3, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        exp1[i] = exp2[i];
    mpn_addmul_1(exp1, exp3, N, scalar);
}

MPOLY_INLINE
void mpoly_monomial_madd_inplace_mp(ulong * exp12, ulong scalar,
                                                   const ulong * exp3, slong N)
{
    mpn_addmul_1(exp12, exp3, N, scalar);
}

MPOLY_INLINE
void mpoly_monomial_msub(ulong * exp1, const ulong * exp2, ulong scalar,
                                                   const ulong * exp3, slong N)
{
   slong i;
   for (i = 0; i < N; i++)
      exp1[i] = exp2[i] - scalar*exp3[i];
}

MPOLY_INLINE
void mpoly_monomial_msub_mp(ulong * exp1, const ulong * exp2, ulong scalar,
                                                   const ulong * exp3, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        exp1[i] = exp2[i];
    mpn_submul_1(exp1, exp3, N, scalar);
}


MPOLY_INLINE
void mpoly_monomial_max(ulong * exp1, const ulong * exp2, const ulong * exp3,
                                               slong bits, slong N, ulong mask)
{
    ulong i, s, m;
    for (i = 0; i < N; i++)
    {
        s = mask + exp2[i] - exp3[i];
        m = mask & s;
        m = m - (m >> (bits - 1));
        exp1[i] = exp3[i] + (s & m);
    }
}

MPOLY_INLINE
void mpoly_monomial_min(ulong * exp1, const ulong * exp2, const ulong * exp3,
                                               slong bits, slong N, ulong mask)
{
    ulong i, s, m;
    for (i = 0; i < N; i++)
    {
        s = mask + exp2[i] - exp3[i];
        m = mask & s;
        m = m - (m >> (bits - 1));
        exp1[i] = exp2[i] - (s & m);
    }
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
int mpoly_monomial_overflows_mp(ulong * exp_ptr, slong N, mp_bitcnt_t bits)
{
    slong i = bits/FLINT_BITS - 1;
    do {
        if ((slong)(exp_ptr[i]) < 0)
            return 1;
        i += bits/FLINT_BITS;
    } while (i < N);

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
int mpoly_monomial_divides_mp(ulong * exp_ptr, const ulong * exp2,
                                 const ulong * exp3, slong N, mp_bitcnt_t bits)
{
    slong i;

    mpn_sub_n(exp_ptr, exp2, exp3, N);

    i = bits/FLINT_BITS - 1;
    do {
        if ((slong)(exp_ptr[i]) < 0)
            return 0;
        i += bits/FLINT_BITS;
    } while (i < N);

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
void mpoly_monomial_mul_ui_mp(ulong * exp2, const ulong * exp3, slong N, ulong c)
{
    mpn_mul_1(exp2, exp3, N, c);
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
                                                slong N, const ulong * cmpmask)
{
    slong i = N - 1;
    do
    {
        if (exp2[i] != exp3[i])
            return (exp3[i]^cmpmask[i]) < (exp2[i]^cmpmask[i]);
    } while (--i >= 0);
    return 0;
}

MPOLY_INLINE
int mpoly_monomial_gt(const ulong * exp2, const ulong * exp3,
                                                slong N, const ulong * cmpmask)
{
    slong i = N - 1;
    do
    {
        if (exp2[i] != exp3[i])
            return (exp3[i]^cmpmask[i]) > (exp2[i]^cmpmask[i]);
    } while (--i >= 0);
    return 0;
}

MPOLY_INLINE
int mpoly_monomial_lt_nomask(const ulong * exp2, const ulong * exp3, slong N)
{
    slong i = N - 1;
    do
    {
        if (exp2[i] != exp3[i])
            return exp2[i] < exp3[i];
    } while (--i >= 0);
    return 0;
}
 
MPOLY_INLINE
int mpoly_monomial_gt_nomask(const ulong * exp2, const ulong * exp3, slong N)
{
    slong i = N - 1;
    do
    {
        if (exp2[i] != exp3[i])
            return exp2[i] > exp3[i];
    } while (--i >= 0);
    return 0;
}

MPOLY_INLINE
int mpoly_monomial_cmp(const ulong * exp2, const ulong * exp3,
                                                slong N, const ulong * cmpmask)
{
    slong i = N - 1;
    do
    {
        if (exp2[i] != exp3[i])
        {
            if ((exp2[i]^cmpmask[i]) > (exp3[i]^cmpmask[i]))
                return 1;
            else
                return -1;
        }
    } while (--i >= 0);
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

FLINT_DLL void mpoly_get_cmpmask(ulong * cmpmask, slong N, slong bits,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_get_ovfmask(ulong * ovfmask, slong N, slong bits,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL mp_bitcnt_t mpoly_exp_bits_required_ui(const ulong * user_exp,
                                                       const mpoly_ctx_t mctx);
FLINT_DLL mp_bitcnt_t mpoly_exp_bits_required_ffmpz(const fmpz * user_exp,
                                                       const mpoly_ctx_t mctx);
FLINT_DLL mp_bitcnt_t mpoly_exp_bits_required_pfmpz(fmpz * const * user_exp,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL slong mpoly_fix_bits(slong bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_pack_vec_ui(ulong * exp1, const ulong * exp2, slong bits,
                                                     slong nfields, slong len);
FLINT_DLL void mpoly_pack_vec_fmpz(ulong * exp1, const fmpz * exp2, mp_bitcnt_t bits,
                                                     slong nfields, slong len);

FLINT_DLL void mpoly_unpack_vec_ui(ulong * exp1, const ulong * exp2, slong bits,
                                                     slong nfields, slong len);
FLINT_DLL void mpoly_unpack_vec_fmpz(fmpz * exp1, const ulong * exp2, mp_bitcnt_t bits,
                                                     slong nfields, slong len);

FLINT_DLL void mpoly_get_monomial_ui(ulong * exps, const ulong * poly_exps,
                                           slong bits, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_get_monomial_ffmpz(fmpz * exps, const ulong * poly_exps,
                                     mp_bitcnt_t bits, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_get_monomial_pfmpz(fmpz ** exps, const ulong * poly_exps,
                                     mp_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_set_monomial_ui(ulong * exp1, const ulong * exp2,
                                           slong bits, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_set_monomial_ffmpz(ulong * exp1, const fmpz * exp2,
                                     mp_bitcnt_t bits, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_set_monomial_pfmpz(ulong * exp1, fmpz * const * exp2,
                                     mp_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_repack_monomials(ulong * exps1, slong bits1,
                                const ulong * exps2, slong bits2, slong len,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_pack_monomials_tight(ulong * exp1,
                  const ulong * exp2, slong len, const slong * mults, 
                                                        slong num, slong bits);

FLINT_DLL void mpoly_unpack_monomials_tight(ulong * e1, ulong * e2, slong len,
                                         slong * mults, slong num, slong bits);

FLINT_DLL int mpoly_monomial_exists(slong * index, const ulong * poly_exps,
                 const ulong * exp, slong len, slong N, const ulong * cmpmask);

FLINT_DLL void mpoly_gen_fields_ui(ulong * exp, slong var, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_gen_fields_fmpz(fmpz * exp, slong var, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_max_fields_ui(ulong * max_fields, const ulong * poly_exps,
                                slong len, slong bits, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_max_fields_fmpz(fmpz * max_fields, const ulong * poly_exps,
                                slong len, slong bits, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_degrees_fit_si(const ulong * poly_exps,
                                slong len, slong bits, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_degrees_si(slong * user_degs, const ulong * poly_exps,
                                slong len, slong bits, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_degrees_ffmpz(fmpz * user_degs, const ulong * poly_exps,
                                slong len, slong bits, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_degrees_pfmpz(fmpz ** user_degs, const ulong * poly_exps,
                                slong len, slong bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_search_monomials(
                slong ** e_ind, ulong * e, slong * e_score,
                slong * t1, slong * t2, slong *t3,
                slong lower, slong upper,
                const ulong * a, slong a_len, const ulong * b, slong b_len,
                                               slong N, const ulong * cmpmask);

FLINT_DLL void mpoly_main_variable_split_LEX(slong * ind, ulong * pexp,
                                 const ulong * Aexp, slong l1, slong Alen,
                                  const ulong * mults, slong num, slong Abits);

FLINT_DLL void mpoly_main_variable_split_DEG(slong * ind, ulong * pexp,
                                  const ulong * Aexp,  slong l1, slong Alen,
                                            ulong deg, slong num, slong Abits);

FLINT_DLL int mpoly_termexp_fits_si(ulong * exps, slong bits,
                                              slong n, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_termexp_fits_ui(ulong * exps, slong bits,
                                              slong n, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_monomials_valid_test(ulong * exps, slong len, slong bits, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_monomials_overflow_test(ulong * exps, slong len, slong bits, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_monomials_inorder_test(ulong * exps, slong len, slong bits, const mpoly_ctx_t mctx);

/* info related to zippel interpolation **************************************/
typedef struct
{
    slong nvars;
    slong * Adegs;
    slong * Bdegs;
    slong * perm;
} mpoly_zipinfo_struct;
typedef mpoly_zipinfo_struct mpoly_zipinfo_t[1];

void mpoly_zipinfo_init(mpoly_zipinfo_t zinfo, slong nvars);

void mpoly_zipinfo_clear(mpoly_zipinfo_t zinfo);

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
                                                         const ulong * cmpmask)
{
   ulong * exp;
   slong i, j, s = --(*heap_len);
   mpoly_heap_t * x = heap[1].next;

   i = 1;
   j = 2;

   while (j < s)
   {
      if (!mpoly_monomial_gt(heap[j + 1].exp, heap[j].exp, N, cmpmask))
         j++;
      heap[i] = heap[j];
      i = j;
      j = HEAP_LEFT(j);
   }

   /* insert last element into heap[i] */
   exp = heap[s].exp;
   j = HEAP_PARENT(i);

   while (i > 1 && mpoly_monomial_gt(heap[j].exp, exp, N, cmpmask))
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
       slong * next_loc, slong * heap_len, slong N, const ulong * cmpmask)
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
      if (!mpoly_monomial_gt(heap[j].exp, exp, N, cmpmask))
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
   slong shift = bits*(k - 1);

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
