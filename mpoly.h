/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017-2019 Daniel Schultz

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

#include "string.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "ulong_extras.h"
#include "thread_pool.h"

#ifdef __cplusplus
 extern "C" {
#endif


#define MPOLY_MIN_BITS (UWORD(8))    /* minimum number of bits to pack into */

/* choose m so that (m + 1)/(n - m) ~= la/lb, i.e. m = (n*la - lb)/(la + lb) */
MPOLY_INLINE slong mpoly_divide_threads(slong n, double la, double lb)
{
    double m_double = (n*la - lb)/(la + lb);
    slong m = m_double + (2*m_double > n ? -0.5 : 0.5);

    /* input must satisfy */
    FLINT_ASSERT(n > 0);

    if (m <= 0)
        m = 0;

    if (m >= n - 1)
        m = n - 1;

    /* output must satisfy */
    FLINT_ASSERT(m >= 0);
    FLINT_ASSERT(m < n);
    return m;
}


typedef enum {
    ORD_LEX,
    ORD_DEGLEX,
    ORD_DEGREVLEX
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
    slong lut_words_per_exp[FLINT_BITS];
    unsigned char lut_fix_bits[FLINT_BITS]; /* FLINT_BITS < 256 */
} mpoly_ctx_struct;

typedef mpoly_ctx_struct mpoly_ctx_t[1];

FLINT_DLL void mpoly_ctx_init(mpoly_ctx_t ctx, slong nvars, const ordering_t ord);

FLINT_DLL void mpoly_ctx_init_rand(mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars);

FLINT_DLL void mpoly_monomial_randbits_fmpz(fmpz * exp, flint_rand_t state, flint_bitcnt_t exp_bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_ctx_clear(mpoly_ctx_t mctx);

/*
    number of words used by an exponent vector packed into "bits" bits:
    we must have either
        (mp) bits > FLINT_BITS and bits % FLINT_BITS == 0, or
        (sp) MPOLY_MIN_BITS <= bits <= FLINT_BITS
*/
MPOLY_INLINE
slong mpoly_words_per_exp_sp(flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    FLINT_ASSERT(0 < bits);
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(mctx->lut_words_per_exp[bits - 1]
                 == (mctx->nfields - 1)/(FLINT_BITS/bits) + 1);
    return mctx->lut_words_per_exp[bits - 1];
}

MPOLY_INLINE
slong mpoly_words_per_exp_mp(flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    FLINT_ASSERT(bits % FLINT_BITS == 0);
    return bits/FLINT_BITS*mctx->nfields;
}

MPOLY_INLINE
slong mpoly_words_per_exp(flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    if (bits <= FLINT_BITS)
        return mpoly_words_per_exp_sp(bits, mctx);
    else
        return mpoly_words_per_exp_mp(bits, mctx);
}

/*
    If "bits" is simply the number of bits needed to pack an exponent vector,
    possibly upgrade it so that it is either
        (mp) a multiple of FLINT_BITS in the mp case, or
        (sp) as big as possible without increasing words_per_exp in the sp case
    The upgrade in (mp) is manditory, while the upgrade in (sp) is simply nice.
*/
MPOLY_INLINE
flint_bitcnt_t mpoly_fix_bits(flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    FLINT_ASSERT(bits > 0);
    if (bits <= FLINT_BITS)
        return mctx->lut_fix_bits[bits - 1];
    else
        return (bits + FLINT_BITS - 1)/FLINT_BITS*FLINT_BITS;
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

FLINT_DLL mpoly_rbnode_struct * mpoly_rbtree_get(int * new_node,
                                         struct mpoly_rbtree *tree, slong rcx);

FLINT_DLL mpoly_rbnode_struct * mpoly_rbtree_get_fmpz(int * new_node,
                                        struct mpoly_rbtree *tree, fmpz_t rcx);

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
void mpoly_monomial_msub_ui_array(ulong * exp1, const ulong * exp2,
                                     const ulong * scalar, slong scalar_limbs,
                                                   const ulong * exp3, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        exp1[i] = exp2[i];
    FLINT_ASSERT(scalar_limbs <= N);
    for (i = 0; i < scalar_limbs; i++)
        mpn_submul_1(exp1 + i, exp3, N - i, scalar[i]);
}

MPOLY_INLINE
void mpoly_monomial_madd_ui_array(ulong * exp1, const ulong * exp2,
                                     const ulong * scalar, slong scalar_limbs,
                                                   const ulong * exp3, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        exp1[i] = exp2[i];
    FLINT_ASSERT(scalar_limbs <= N);
    for (i = 0; i < scalar_limbs; i++)
        mpn_addmul_1(exp1 + i, exp3, N - i, scalar[i]);
}

MPOLY_INLINE
void mpoly_monomial_madd_fmpz(ulong * exp1, const ulong * exp2,
                             const fmpz_t scalar, const ulong * exp3, slong N)
{
    if (COEFF_IS_MPZ(*scalar))
    {
        __mpz_struct * mpz = COEFF_TO_PTR(*scalar);
        mpoly_monomial_madd_ui_array(exp1, exp2,
                                           mpz->_mp_d, mpz->_mp_size, exp3, N);
    }
    else
    {
        mpoly_monomial_madd_mp(exp1, exp2, *scalar, exp3, N);
    }
}

MPOLY_INLINE
ulong mpoly_overflow_mask_sp(flint_bitcnt_t bits)
{
    slong i;
    ulong mask = 0;

    FLINT_ASSERT(bits <= FLINT_BITS);

    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    return mask;
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
void mpoly_monomial_max_mp(ulong * exp1, const ulong * exp2, const ulong * exp3,
                                                     flint_bitcnt_t bits, slong N)
{
    slong i, j;
    for (i = 0; i < N; i += bits/FLINT_BITS)
    {
        const ulong * t = exp2;
        for (j = bits/FLINT_BITS - 1; j >= 0; j--)
        {
            if (exp3[i + j] != exp2[i + j])
            {
                if (exp3[i + j] > exp2[i + j])
                    t = exp3;
                break;
            }
        }
        for (j = 0; j < bits/FLINT_BITS; j++)
        {
            exp1[i + j] = t[i + j];
        }
    }
}

MPOLY_INLINE
void mpoly_monomial_min_mp(ulong * exp1, const ulong * exp2, const ulong * exp3,
                                                     flint_bitcnt_t bits, slong N)
{
    slong i, j;
    for (i = 0; i < N; i += bits/FLINT_BITS)
    {
        const ulong * t = exp2;
        for (j = bits/FLINT_BITS - 1; j >= 0; j--)
        {
            if (exp3[i + j] != exp2[i + j])
            {
                if (exp3[i + j] < exp2[i + j])
                    t = exp3;
                break;
            }
        }
        for (j = 0; j < bits/FLINT_BITS; j++)
        {
            exp1[i + j] = t[i + j];
        }
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
int mpoly_monomial_overflows_mp(ulong * exp_ptr, slong N, flint_bitcnt_t bits)
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
                                 const ulong * exp3, slong N, flint_bitcnt_t bits)
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
int mpoly_monomial_divides_mp_test(const ulong * exp2,
                                 const ulong * exp3, slong N, flint_bitcnt_t bits)
{
    slong i, j;

    i = 0;
    do {
        for (j = bits/FLINT_BITS - 1; j >= 0; j--)
        {
            if (exp2[i + j] > exp3[i + j])
                break;
            if (exp2[i + j] < exp3[i + j])
                return 0;
        }
        i += bits/FLINT_BITS;
    } while (i < N);

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
void mpoly_monomial_set_extra(ulong * exp2, const ulong * exp3,
                                            slong N, slong offset, ulong extra)
{
    slong i;
    for (i = 0; i < N; i++)
    {
        exp2[i] = exp3[i] + (i == offset ? extra : 0);
    }
}

MPOLY_INLINE
void mpoly_copy_monomials(ulong * exp1, const ulong * exp2, slong len, slong N)
{
    memcpy(exp1, exp2, N*len*sizeof(ulong));
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
void mpoly_monomial_mul_ui(ulong * exp2, const ulong * exp3, slong N, ulong c)
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

FLINT_DLL void mpoly_monomial_mul_fmpz(ulong * exp2, const ulong * exp3,
                                                            slong N, fmpz_t c);

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
int mpoly_monomial_equal_extra(const ulong * exp2, const ulong * exp3,
                                            slong N, slong offset, ulong extra)
{
   slong i;

   for (i = 0; i < N; i++)
   {
      ulong e3 = exp3[i] + ((i == offset) ? extra : 0);
      if (exp2[i] != e3)
         return 0;
   }

   return 1;
}

MPOLY_INLINE
int mpoly_monomial_cmp1(ulong a, ulong b, ulong cmpmask)
{
    if ((a^cmpmask) != (b^cmpmask))
    {
        if ((a^cmpmask) > (b^cmpmask))
            return 1;
        else
            return -1;
    }
    return 0;
}

MPOLY_INLINE
int mpoly_monomial_gt1(ulong a, ulong b, ulong cmpmask)
{
    return (a^cmpmask) > (b^cmpmask);
}

MPOLY_INLINE
int mpoly_monomial_ge1(ulong a, ulong b, ulong cmpmask)
{
    return (a^cmpmask) >= (b^cmpmask);
}

MPOLY_INLINE
int mpoly_monomial_lt(const ulong * exp3, const ulong * exp2,
                                                slong N, const ulong * cmpmask)
{
    slong i = N - 1;
    do {
        if (exp2[i] != exp3[i])
        {
            return (exp3[i]^cmpmask[i]) < (exp2[i]^cmpmask[i]);
        }
    } while (--i >= 0);
    return 0;
}

MPOLY_INLINE
int mpoly_monomial_gt(const ulong * exp3, const ulong * exp2,
                                                slong N, const ulong * cmpmask)
{
    slong i = N - 1;
    do {
        if (exp2[i] != exp3[i])
        {
            return (exp3[i]^cmpmask[i]) > (exp2[i]^cmpmask[i]);
        }
    } while (--i >= 0);
    return 0;
}

MPOLY_INLINE
int mpoly_monomial_lt_nomask(const ulong * exp2, const ulong * exp3, slong N)
{
    slong i = N - 1;
    do {
        if (exp2[i] != exp3[i])
        {
            return exp2[i] < exp3[i];
        }
    } while (--i >= 0);
    return 0;
}
 
MPOLY_INLINE
int mpoly_monomial_gt_nomask(const ulong * exp2, const ulong * exp3, slong N)
{
    slong i = N - 1;
    do {
        if (exp2[i] != exp3[i])
        {
            return exp2[i] > exp3[i];
        }
    } while (--i >= 0);
    return 0;
}

MPOLY_INLINE
int mpoly_monomial_lt_nomask_extra(const ulong * exp2, const ulong * exp3,
                                            slong N, slong offset, ulong extra)
{
    slong i = N - 1;
    do {
        ulong e3 = exp3[i] + ((i == offset) ? extra : 0);
        if (exp2[i] != e3)
        {
            return exp2[i] < e3;
        }
    } while (--i >= 0);
    return 0;
}
 
MPOLY_INLINE
int mpoly_monomial_gt_nomask_extra(const ulong * exp2, const ulong * exp3,
                                            slong N, slong offset, ulong extra)
{
    slong i = N - 1;
    do {
        ulong e3 = exp3[i] + ((i == offset) ? extra : 0);
        if (exp2[i] != e3)
        {
            return exp2[i] > e3;
        }
    } while (--i >= 0);
    return 0;
}

MPOLY_INLINE
int mpoly_monomial_cmp(const ulong * exp2, const ulong * exp3,
                                                slong N, const ulong * cmpmask)
{
    slong i = N - 1;
    do {
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
int mpoly_monomial_cmp_nomask(const ulong * exp2, const ulong * exp3, slong N)
{
    slong i = N - 1;
    do {
        if (exp2[i] != exp3[i])
        {
            if (exp2[i] > exp3[i])
                return 1;
            else
                return -1;
        }
    } while (--i >= 0);
    return 0;
}

MPOLY_INLINE
int mpoly_monomial_cmp_nomask_extra(const ulong * exp2, const ulong * exp3,
                                            slong N, slong offset, ulong extra)
{
    slong i = N - 1;
    do {
        ulong e3 = exp3[i] + ((i == offset) ? extra : 0);
        if (exp2[i] != e3)
        {
            if (exp2[i] > e3)
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


/* generators ****************************************************************/

FLINT_DLL void mpoly_gen_fields_ui(ulong * exp, slong var,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_gen_fields_fmpz(fmpz * exp, slong var,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL flint_bitcnt_t mpoly_gen_bits_required(slong var, const mpoly_ctx_t mctx);

/* return the index in the fields where the generator of index v is stored */
MPOLY_INLINE slong mpoly_gen_index(slong v, const mpoly_ctx_t mctx)
{
    return mctx->rev ? v : mctx->nvars - 1 - v;
}

FLINT_DLL void mpoly_gen_offset_shift_sp(slong * offset, slong * shift,
                          slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_gen_monomial_offset_shift_sp(ulong * mexp, slong * offset,
           slong * shift, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_gen_monomial_sp(ulong * oneexp, slong var,
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL slong mpoly_gen_offset_mp(slong var, 
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL slong mpoly_gen_monomial_offset_mp(ulong * mexp, slong var,
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void fmpz_mat_mul_vec(fmpz * v, const fmpz_mat_t M, fmpz * u);

FLINT_DLL void mpoly_compose_mat_gen(fmpz_mat_t M, const slong * c,
                            const mpoly_ctx_t mctxB, const mpoly_ctx_t mctxAC);

FLINT_DLL void mpoly_compose_mat_fill_column(fmpz_mat_t M, const ulong * Cexp,
                    flint_bitcnt_t Cbits, slong Bvar, const mpoly_ctx_t mctxB,
                                                     const mpoly_ctx_t mctxAC);

/* Monomial arrays ***********************************************************/

FLINT_DLL void mpoly_get_cmpmask(ulong * cmpmask, slong N, flint_bitcnt_t bits,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_get_ovfmask(ulong * ovfmask, slong N, flint_bitcnt_t bits,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL flint_bitcnt_t mpoly_exp_bits_required_ui(const ulong * user_exp,
                                                       const mpoly_ctx_t mctx);
FLINT_DLL flint_bitcnt_t mpoly_exp_bits_required_ffmpz(const fmpz * user_exp,
                                                       const mpoly_ctx_t mctx);
FLINT_DLL flint_bitcnt_t mpoly_exp_bits_required_pfmpz(fmpz * const * user_exp,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_pack_vec_ui(ulong * exp1, const ulong * exp2,
                                flint_bitcnt_t bits, slong nfields, slong len);
FLINT_DLL void mpoly_pack_vec_fmpz(ulong * exp1, const fmpz * exp2,
                                flint_bitcnt_t bits, slong nfields, slong len);

FLINT_DLL void mpoly_unpack_vec_ui(ulong * exp1, const ulong * exp2,
                                flint_bitcnt_t bits, slong nfields, slong len);

FLINT_DLL void mpoly_unpack_vec_fmpz(fmpz * exp1, const ulong * exp2,
                                flint_bitcnt_t bits, slong nfields, slong len);

FLINT_DLL void mpoly_get_monomial_ui_unpacked_ffmpz(ulong * user_exps,
                               const fmpz * poly_exps, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_get_monomial_ffmpz_unpacked_ffmpz(fmpz * user_exps,
                               const fmpz * poly_exps, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_get_monomial_pfmpz_unpacked_ffmpz(fmpz ** user_exps,
                               const fmpz * poly_exps, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_get_monomial_ui_unpacked_ui(ulong * user_exps,
                              const ulong * poly_exps, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_get_monomial_ui_sp(ulong * user_exps,
         const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_get_monomial_ui_mp(ulong * user_exps,
         const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_get_monomial_si_mp(slong * user_exps,
         const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

MPOLY_INLINE void mpoly_get_monomial_ui(ulong * user_exps,
          const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    if (bits <= FLINT_BITS)
        mpoly_get_monomial_ui_sp(user_exps, poly_exps, bits, mctx);
    else
        mpoly_get_monomial_ui_mp(user_exps, poly_exps, bits, mctx);
}

MPOLY_INLINE void mpoly_get_monomial_si(slong * user_exps,
          const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
/* if bits <= FLINT_BITS and poly_exps is canonical, everything should be ok */
    if (bits <= FLINT_BITS)
        mpoly_get_monomial_ui_sp((ulong *) user_exps, poly_exps, bits, mctx);
    else
        mpoly_get_monomial_si_mp(user_exps, poly_exps, bits, mctx);
}

FLINT_DLL ulong mpoly_get_monomial_var_exp_ui_sp(const ulong * poly_exps,
                       slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL ulong mpoly_get_monomial_var_exp_ui_mp(const ulong * poly_exps,
                       slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL slong mpoly_get_monomial_var_exp_si_mp(const ulong * poly_exps,
                       slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

MPOLY_INLINE ulong mpoly_get_monomial_var_exp_ui(const ulong * poly_exps,
                        slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    if (bits <= FLINT_BITS)
        return mpoly_get_monomial_var_exp_ui_sp(poly_exps, var, bits, mctx);
    else
        return mpoly_get_monomial_var_exp_ui_mp(poly_exps, var, bits, mctx);
}

MPOLY_INLINE slong mpoly_get_monomial_var_exp_si(const ulong * poly_exps,
                        slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    if (bits <= FLINT_BITS)
        return (slong) mpoly_get_monomial_var_exp_ui_sp(poly_exps, var, bits, mctx);
    else
        return mpoly_get_monomial_var_exp_si_mp(poly_exps, var, bits, mctx);
}

FLINT_DLL void mpoly_get_monomial_ffmpz(fmpz * exps, const ulong * poly_exps,
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_get_monomial_pfmpz(fmpz ** exps, const ulong * poly_exps,
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_set_monomial_ui(ulong * exp1, const ulong * exp2,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_set_monomial_ffmpz(ulong * exp1, const fmpz * exp2,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_set_monomial_pfmpz(ulong * exp1, fmpz * const * exp2,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_repack_monomials(ulong * exps1, flint_bitcnt_t bits1,
                        const ulong * exps2, flint_bitcnt_t bits2, slong len,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_pack_monomials_tight(ulong * exp1,
                  const ulong * exp2, slong len, const slong * mults, 
                                                        slong num, slong bits);

FLINT_DLL void mpoly_unpack_monomials_tight(ulong * e1, ulong * e2, slong len,
                                         slong * mults, slong num, slong bits);

FLINT_DLL int mpoly_monomial_exists(slong * index, const ulong * poly_exps,
                 const ulong * exp, slong len, slong N, const ulong * cmpmask);

FLINT_DLL slong mpoly_monomial_index_ui(const ulong * Aexp, flint_bitcnt_t Abits,
                     slong Alength, const ulong * exp, const mpoly_ctx_t mctx);

FLINT_DLL slong mpoly_monomial_index_pfmpz(const ulong * Aexp, flint_bitcnt_t Abits,
                    slong Alength, fmpz * const * exp, const mpoly_ctx_t mctx);

FLINT_DLL slong mpoly_monomial_index_monomial(const ulong * Aexp,
                      flint_bitcnt_t Abits, slong Alength, const ulong * Mexp,
                                    flint_bitcnt_t Mbits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_min_fields_ui_sp(ulong * min_fields, const ulong * poly_exps,
                          slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_min_fields_fmpz(fmpz * min_fields, const ulong * poly_exps,
                          slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_max_fields_ui_sp(ulong * max_fields, const ulong * poly_exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_max_fields_fmpz(fmpz * max_fields, const ulong * poly_exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_degrees_fit_si(const ulong * poly_exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_degrees_si(slong * user_degs, const ulong * poly_exps,
                          slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_degrees_si_threaded(slong * user_degs, const ulong * poly_exps,
                         slong len,  flint_bitcnt_t bits, const mpoly_ctx_t mctx,
                        const thread_pool_handle * handles, slong num_handles);

FLINT_DLL void mpoly_degrees_ffmpz(fmpz * user_degs, const ulong * poly_exps,
                          slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_degrees_pfmpz(fmpz ** user_degs, const ulong * poly_exps,
                          slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL slong mpoly_degree_si(const ulong * poly_exps,
               slong len, flint_bitcnt_t bits, slong var, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_degree_fmpz(fmpz_t deg, const ulong * poly_exps,
               slong len, flint_bitcnt_t bits, slong var, const mpoly_ctx_t mctx);

FLINT_DLL int  mpoly_total_degree_fits_si(const ulong * exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL slong mpoly_total_degree_si(const ulong * exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_total_degree_fmpz(fmpz_t totdeg, const ulong * exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_total_degree_fmpz_ref(fmpz_t totdeg, const ulong * exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

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

FLINT_DLL int mpoly_term_exp_fits_si(ulong * exps, flint_bitcnt_t bits,
                                              slong n, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_term_exp_fits_ui(ulong * exps, flint_bitcnt_t bits,
                                              slong n, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_is_gen(ulong * exps, slong var,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_monomials_valid_test(ulong * exps, slong len,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_monomials_overflow_test(ulong * exps, slong len,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_monomials_inorder_test(ulong * exps, slong len,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_reverse(ulong * Aexp, const ulong * Bexp, slong len, slong N);

FLINT_DLL void mpoly_monomials_deflation(fmpz * shift, fmpz * stride,
                        const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_monomials_deflate(ulong * Aexps, flint_bitcnt_t Abits,
                        const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength,
              const fmpz * shift, const fmpz * stride, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_monomials_inflate(ulong * Aexps, flint_bitcnt_t Abits,
                        const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength,
              const fmpz * shift, const fmpz * stride, const mpoly_ctx_t mctx);

FLINT_DLL void _mpoly_gen_shift_right(ulong * Aexp, flint_bitcnt_t Abits,
               slong Alength, slong var, ulong amount, const mpoly_ctx_t mctx);

FLINT_DLL void _mpoly_gen_shift_left(ulong * Aexp, flint_bitcnt_t Abits,
               slong Alength, slong var, ulong amount, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_monomial_cmp_general(ulong * Aexp, flint_bitcnt_t Abits,
                      ulong * Bexp, flint_bitcnt_t Bbits, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_monomials_shift_right_ui(ulong * Aexps, flint_bitcnt_t Abits,
               slong Alength, const ulong * user_exps, const mpoly_ctx_t mctx);

FLINT_DLL int mpoly_monomial_cofactors(fmpz * Abarexps, fmpz * Bbarexps,
                                    const ulong * Aexps, flint_bitcnt_t Abits,
                                    const ulong * Bexps, flint_bitcnt_t Bbits,
                                        slong length,  const mpoly_ctx_t mctx);

/* info related to gcd calculation *******************************************/

typedef struct
{
    ulong * Amax_exp;
    ulong * Amin_exp;
    ulong * Astride;
    slong * Adeflate_deg;
    slong * Alead_count;
    slong * Atail_count;

    ulong * Bmax_exp;
    ulong * Bmin_exp;
    ulong * Bstride;
    slong * Bdeflate_deg;
    slong * Blead_count;
    slong * Btail_count;

    ulong * Gmin_exp;
    ulong * Abarmin_exp;
    ulong * Bbarmin_exp;
    ulong * Gstride;
    slong * Gterm_count_est;
    slong * Gdeflate_deg_bound;

    flint_bitcnt_t Gbits, Abarbits, Bbarbits;

    slong mvars;

    double Adensity;
    double Bdensity;

    double brown_time_est, bma_time_est, zippel_time_est;
    slong * brown_perm, * bma_perm, * zippel_perm;
    int can_use_brown, can_use_bma, can_use_zippel;
    int Gdeflate_deg_bounds_are_nice; /* all of Gdeflate_deg_bound came from real gcd computations */

    char * data;
} mpoly_gcd_info_struct;

typedef mpoly_gcd_info_struct mpoly_gcd_info_t[1];

FLINT_DLL void mpoly_gcd_info_init(mpoly_gcd_info_t I, slong nvars);

FLINT_DLL void mpoly_gcd_info_clear(mpoly_gcd_info_t I);

FLINT_DLL void mpoly_gcd_info_limits(ulong * Amax_exp, ulong * Amin_exp,
                       slong * Amax_exp_count, slong * Amin_exp_count,
                       const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                                                       const mpoly_ctx_t mctx);
FLINT_DLL void mpoly_gcd_info_stride(ulong * strides,
          const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                             const ulong * Amax_exp, const ulong * Amin_exp,
          const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength,
                             const ulong * Bmax_exp, const ulong * Bmin_exp,
                                                       const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_gcd_info_set_perm(mpoly_gcd_info_t I,
                         slong Alength, slong Blength, const mpoly_ctx_t mctx);

FLINT_DLL slong mpoly_gcd_info_get_brown_upper_limit(const mpoly_gcd_info_t I,
                                                       slong var, slong bound);

FLINT_DLL void mpoly_gcd_info_measure_brown(mpoly_gcd_info_t I,
                         slong Alength, slong Blength, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_gcd_info_measure_bma(mpoly_gcd_info_t I,
                         slong Alength, slong Blength, const mpoly_ctx_t mctx);

FLINT_DLL void mpoly_gcd_info_measure_zippel(mpoly_gcd_info_t I,
                         slong Alength, slong Blength, const mpoly_ctx_t mctx);

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

void _fmpz_vec_content_chained(fmpz_t res, const fmpz * vec, slong len);

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
      ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[1].next;
      heap[1].next = x;

      return;
   }

   if (*next_loc < *heap_len)
   {
      if (exp == heap[*next_loc].exp)
      {
         ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[*next_loc].next;
         heap[*next_loc].next = x;
         return;
      }
   }



   while ((j = HEAP_PARENT(i)) >= 1)
   {
      if (exp == heap[j].exp)
      {
         ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[j].next;
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
   mpoly_heap_t * x = (mpoly_heap_t *) heap[1].next;

   i = 1;
   j = 2;

   while (j < s)
   {
      if (!mpoly_monomial_gt(heap[j].exp, heap[j + 1].exp, N, cmpmask))
         j++;
      heap[i] = heap[j];
      i = j;
      j = HEAP_LEFT(j);
   }

   /* insert last element into heap[i] */
   exp = heap[s].exp;
   j = HEAP_PARENT(i);

   while (i > 1 && mpoly_monomial_gt(exp, heap[j].exp, N, cmpmask))
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
      ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[1].next;
      heap[1].next = x;

      return 0;
   }

   if (*next_loc < *heap_len)
   {
      if (mpoly_monomial_equal(exp, heap[*next_loc].exp, N))
      {
         ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[*next_loc].next;
         heap[*next_loc].next = x;
         return 0;
      }
   }

   while ((j = HEAP_PARENT(i)) >= 1)
   {
      if (!mpoly_monomial_gt(exp, heap[j].exp, N, cmpmask))
         break;

      i = j;
   }

   if (j >= 1 && mpoly_monomial_equal(exp, heap[j].exp, N))
   {
      ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[j].next;
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
