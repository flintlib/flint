/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPOLY_H
#define MPOLY_H

#ifdef MPOLY_INLINES_C
#define MPOLY_INLINE
#else
#define MPOLY_INLINE static inline
#endif

#include <string.h>
#include "thread_pool.h"
#include "mpoly_types.h"
#include "fmpz.h"

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


/* context *******************************************************************/

void mpoly_ctx_init(mpoly_ctx_t ctx, slong nvars, const ordering_t ord);

void mpoly_ctx_init_rand(mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars);

void mpoly_monomial_randbits_fmpz(fmpz * exp, flint_rand_t state, flint_bitcnt_t exp_bits, const mpoly_ctx_t mctx);

void mpoly_ctx_clear(mpoly_ctx_t mctx);

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

/* red-black with ui keys */
typedef struct {
    ulong key;
    slong up;
    slong left;
    slong right;
    int color;
} mpoly_rbnode_ui_struct;

typedef struct {
    slong length;
    mpoly_rbnode_ui_struct * nodes;
    slong node_alloc;
    char * data;
    slong data_alloc;
    slong data_size;
} mpoly_rbtree_ui_struct;

typedef mpoly_rbtree_ui_struct mpoly_rbtree_ui_t[1];

void mpoly_rbtree_ui_init(mpoly_rbtree_ui_t T, slong data_size);

void mpoly_rbtree_ui_clear(mpoly_rbtree_ui_t T);

void * mpoly_rbtree_ui_lookup(mpoly_rbtree_ui_t T, int * its_new,
                                                                    ulong key);

MPOLY_INLINE slong mpoly_rbtree_ui_head(const mpoly_rbtree_ui_t T)
{
    FLINT_ASSERT(T->nodes[1].left >= 0 || T->length < 1);
    return T->nodes[1].left;
}

/* red-black with fmpz keys */
typedef struct {
    fmpz_t key;
    slong up;
    slong left;
    slong right;
    int color;
} mpoly_rbnode_fmpz_struct;

typedef struct {
    slong length;
    mpoly_rbnode_fmpz_struct * nodes;
    slong node_alloc;
    char * data;
    slong data_alloc;
    slong data_size;
} mpoly_rbtree_fmpz_struct;

typedef mpoly_rbtree_fmpz_struct mpoly_rbtree_fmpz_t[1];

void mpoly_rbtree_fmpz_init(mpoly_rbtree_fmpz_t T, slong data_size);

void mpoly_rbtree_fmpz_clear(mpoly_rbtree_fmpz_t T);

void * mpoly_rbtree_fmpz_lookup(mpoly_rbtree_fmpz_t T, int * its_new,
                                                             const fmpz_t key);

MPOLY_INLINE slong mpoly_rbtree_fmpz_head(const mpoly_rbtree_fmpz_t T)
{
    FLINT_ASSERT(T->nodes[1].left >= 0 || T->length < 1);
    return T->nodes[1].left;
}

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

void mpoly_ordering_print(ordering_t ord);

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
    FLINT_ASSERT(N > 0);
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
    {
        FLINT_ASSERT(N > i);
        mpn_submul_1(exp1 + i, exp3, N - i, scalar[i]);
    }
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

/* mask with high bit set in each field of exponent vector */
MPOLY_INLINE
ulong mpoly_overflow_mask_sp(flint_bitcnt_t bits)
{
    flint_bitcnt_t i;
    ulong mask = 0;

    FLINT_ASSERT(bits <= FLINT_BITS);

    mask = (UWORD(1) << (bits - 1));
    for (i = bits; i < FLINT_BITS; i += bits)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    return mask;
}

MPOLY_INLINE
ulong mpoly_monomial_max1(ulong exp2, ulong exp3,
                                               flint_bitcnt_t bits, ulong mask)
{
    ulong s, m, exp1;
    s = mask + exp2 - exp3;
    m = mask & s;
    m = m - (m >> (bits - 1));
    exp1 = exp3 + (s & m);
    return exp1;
}

MPOLY_INLINE
void mpoly_monomial_max(ulong * exp1, const ulong * exp2, const ulong * exp3,
                                      flint_bitcnt_t bits, slong N, ulong mask)
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
ulong mpoly_monomial_min1(ulong exp2, ulong exp3,
                                               flint_bitcnt_t bits, ulong mask)
{
    ulong s, m, exp1;
    s = mask + exp2 - exp3;
    m = mask & s;
    m = m - (m >> (bits - 1));
    exp1 = exp2 - (s & m);
    return exp1;
}

MPOLY_INLINE
void mpoly_monomial_min(ulong * exp1, const ulong * exp2, const ulong * exp3,
                                      flint_bitcnt_t bits, slong N, ulong mask)
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
int mpoly_monomial_halves(ulong * exp_ptr, const ulong * exp2,
		                                           slong N, ulong mask)
{
    slong i;
    for (i = 0; i < N; i++)
    {
        if (exp2[i] & 1)
            return 0;
        exp_ptr[i] = exp2[i] >> 1;
        if (exp_ptr[i] & mask)
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
int mpoly_monomial_halves_mp(ulong * exp_ptr, const ulong * exp2,
		                                  slong N, flint_bitcnt_t bits)
{
   slong i;
   ulong bw;

   bw = mpn_rshift(exp_ptr, exp2, N, 1);

   if (bw != 0)
      return 0;

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
int mpoly_monomial_halves1(ulong * exp_ptr, const ulong exp2, ulong mask)
{
   if (exp2 & 1)
      return 0;

   (*exp_ptr) = exp2 >> 1;

   if (((exp2 >> 1) & mask) != 0)
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
    if (len > 0)
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
    FLINT_ASSERT(N > 0);
    mpn_mul_1(exp2, exp3, N, c);
}

void mpoly_monomial_mul_fmpz(ulong * exp2, const ulong * exp3,
                                                      slong N, const fmpz_t c);

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


/* ceiling(log_4(x)) - 1 */
MPOLY_INLINE slong mpoly_geobucket_clog4(slong x)
{
    if (x <= 4)
        return 0;
    /*
        FLINT_BIT_COUNT returns unsigned int.
        Signed division is not defined.
        Do the calculation with unsigned ints and then convert to slong.
    */
    return (slong)((FLINT_BIT_COUNT(x - 1) - UWORD(1))/(UWORD(2)));
}

/* single-limb packings ******************************************************/

MPOLY_INLINE
ulong pack_exp2(ulong e0, ulong e1)
{
    return (e0 << (1*(FLINT_BITS/2))) +
           (e1 << (0*(FLINT_BITS/2)));
}

MPOLY_INLINE
ulong pack_exp3(ulong e0, ulong e1, ulong e2)
{
    return (e0 << (2*(FLINT_BITS/3))) +
           (e1 << (1*(FLINT_BITS/3))) +
           (e2 << (0*(FLINT_BITS/3)));
}

MPOLY_INLINE
ulong extract_exp(ulong e, int idx, int nvars)
{
    return (e >> (idx*(FLINT_BITS/nvars))) &
            ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/nvars));
}

ulong _mpoly_bidegree(const ulong * Aexps, flint_bitcnt_t Abits,
                                                       const mpoly_ctx_t mctx);

/* generators ****************************************************************/

void mpoly_gen_fields_ui(ulong * exp, slong var,
                                                       const mpoly_ctx_t mctx);

void mpoly_gen_fields_fmpz(fmpz * exp, slong var,
                                                       const mpoly_ctx_t mctx);

flint_bitcnt_t mpoly_gen_bits_required(slong var, const mpoly_ctx_t mctx);

/* return the index in the fields where the generator of index v is stored */
MPOLY_INLINE slong mpoly_gen_index(slong v, const mpoly_ctx_t mctx)
{
    return mctx->rev ? v : mctx->nvars - 1 - v;
}

void mpoly_gen_offset_shift_sp(slong * offset, slong * shift,
                          slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_gen_monomial_offset_shift_sp(ulong * mexp, slong * offset,
           slong * shift, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_gen_monomial_sp(ulong * oneexp, slong var,
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

slong mpoly_gen_offset_mp(slong var,
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

slong mpoly_gen_monomial_offset_mp(ulong * mexp, slong var,
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void fmpz_mat_mul_vec(fmpz * v, const fmpz_mat_t M, fmpz * u);

void mpoly_compose_mat_gen(fmpz_mat_t M, const slong * c,
                            const mpoly_ctx_t mctxB, const mpoly_ctx_t mctxAC);

void mpoly_compose_mat_fill_column(fmpz_mat_t M, const ulong * Cexp,
                    flint_bitcnt_t Cbits, slong Bvar, const mpoly_ctx_t mctxB,
                                                     const mpoly_ctx_t mctxAC);

/* Monomial arrays ***********************************************************/

void mpoly_get_cmpmask(ulong * cmpmask, slong N, flint_bitcnt_t bits,
                                                       const mpoly_ctx_t mctx);

void mpoly_get_ovfmask(ulong * ovfmask, slong N, flint_bitcnt_t bits,
                                                       const mpoly_ctx_t mctx);

int mpoly_monomials_cmp(const ulong * Aexps, flint_bitcnt_t Abits,
                                  const ulong * Bexps, flint_bitcnt_t Bbits,
                                         slong length, const mpoly_ctx_t mctx);

flint_bitcnt_t mpoly_exp_bits_required_ui(const ulong * user_exp,
                                                       const mpoly_ctx_t mctx);
flint_bitcnt_t mpoly_exp_bits_required_ffmpz(const fmpz * user_exp,
                                                       const mpoly_ctx_t mctx);
flint_bitcnt_t mpoly_exp_bits_required_pfmpz(fmpz * const * user_exp,
                                                       const mpoly_ctx_t mctx);

MPOLY_INLINE
flint_bitcnt_t mpoly_gen_pow_exp_bits_required(slong v, ulong e,
                                                        const mpoly_ctx_t mctx)
{
    return 1 + FLINT_BIT_COUNT(e); /* only lex and deg supported */
}

int mpoly_is_poly(const ulong * Aexps, slong Alen,
                      flint_bitcnt_t Abits, slong var, const mpoly_ctx_t mctx);

void mpoly_pack_vec_ui(ulong * exp1, const ulong * exp2,
                                flint_bitcnt_t bits, slong nfields, slong len);
void mpoly_pack_vec_fmpz(ulong * exp1, const fmpz * exp2,
                                flint_bitcnt_t bits, slong nfields, slong len);

void mpoly_unpack_vec_ui(ulong * exp1, const ulong * exp2,
                                flint_bitcnt_t bits, slong nfields, slong len);

void mpoly_unpack_vec_fmpz(fmpz * exp1, const ulong * exp2,
                                flint_bitcnt_t bits, slong nfields, slong len);

void mpoly_get_monomial_ui_unpacked_ffmpz(ulong * user_exps,
                               const fmpz * poly_exps, const mpoly_ctx_t mctx);

void mpoly_get_monomial_ffmpz_unpacked_ffmpz(fmpz * user_exps,
                               const fmpz * poly_exps, const mpoly_ctx_t mctx);

void mpoly_get_monomial_pfmpz_unpacked_ffmpz(fmpz ** user_exps,
                               const fmpz * poly_exps, const mpoly_ctx_t mctx);

void mpoly_get_monomial_ui_unpacked_ui(ulong * user_exps,
                              const ulong * poly_exps, const mpoly_ctx_t mctx);

void mpoly_get_monomial_ui_sp(ulong * user_exps,
         const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_get_monomial_ui_mp(ulong * user_exps,
         const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_get_monomial_si_mp(slong * user_exps,
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

ulong mpoly_get_monomial_var_exp_ui_sp(const ulong * poly_exps,
                       slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

ulong mpoly_get_monomial_var_exp_ui_mp(const ulong * poly_exps,
                       slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

slong mpoly_get_monomial_var_exp_si_mp(const ulong * poly_exps,
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

void mpoly_get_monomial_ffmpz(fmpz * exps, const ulong * poly_exps,
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_get_monomial_pfmpz(fmpz ** exps, const ulong * poly_exps,
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_set_monomial_ui(ulong * exp1, const ulong * exp2,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);
void mpoly_set_monomial_ffmpz(ulong * exp1, const fmpz * exp2,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);
void mpoly_set_monomial_pfmpz(ulong * exp1, fmpz * const * exp2,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);

int mpoly_repack_monomials(ulong * exps1, flint_bitcnt_t bits1,
                        const ulong * exps2, flint_bitcnt_t bits2, slong len,
                                                       const mpoly_ctx_t mctx);

void mpoly_pack_monomials_tight(ulong * exp1,
                  const ulong * exp2, slong len, const slong * mults,
                                                        slong num, slong bits);

void mpoly_unpack_monomials_tight(ulong * e1, ulong * e2, slong len,
                                         slong * mults, slong num, slong bits);

int mpoly_monomial_exists(slong * index, const ulong * poly_exps,
                 const ulong * exp, slong len, slong N, const ulong * cmpmask);

slong mpoly_monomial_index1_nomask(ulong * Aexps, slong Alen, ulong e);

slong mpoly_monomial_index_ui(const ulong * Aexp, flint_bitcnt_t Abits,
                     slong Alength, const ulong * exp, const mpoly_ctx_t mctx);

slong mpoly_monomial_index_pfmpz(const ulong * Aexp, flint_bitcnt_t Abits,
                    slong Alength, fmpz * const * exp, const mpoly_ctx_t mctx);

slong mpoly_monomial_index_monomial(const ulong * Aexp,
                      flint_bitcnt_t Abits, slong Alength, const ulong * Mexp,
                                    flint_bitcnt_t Mbits, const mpoly_ctx_t mctx);

void mpoly_min_fields_ui_sp(ulong * min_fields, const ulong * poly_exps,
                          slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_min_fields_fmpz(fmpz * min_fields, const ulong * poly_exps,
                          slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_max_fields_ui_sp(ulong * max_fields, const ulong * poly_exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_max_fields_fmpz(fmpz * max_fields, const ulong * poly_exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

int mpoly_degrees_fit_si(const ulong * poly_exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_degrees_si(slong * user_degs, const ulong * poly_exps,
                          slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_degrees_si_threaded(slong * user_degs, const ulong * poly_exps,
                         slong len,  flint_bitcnt_t bits, const mpoly_ctx_t mctx,
                        const thread_pool_handle * handles, slong num_handles);

void mpoly_degrees_ffmpz(fmpz * user_degs, const ulong * poly_exps,
                          slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_degrees_pfmpz(fmpz ** user_degs, const ulong * poly_exps,
                          slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

slong mpoly_degree_si(const ulong * poly_exps,
               slong len, flint_bitcnt_t bits, slong var, const mpoly_ctx_t mctx);

void mpoly_degree_fmpz(fmpz_t deg, const ulong * poly_exps,
               slong len, flint_bitcnt_t bits, slong var, const mpoly_ctx_t mctx);

int  mpoly_total_degree_fits_si(const ulong * exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

slong mpoly_total_degree_si(const ulong * exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_total_degree_fmpz(fmpz_t totdeg, const ulong * exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_total_degree_fmpz_ref(fmpz_t totdeg, const ulong * exps,
                                slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_used_vars_or(int * used, const ulong * exps,
                       slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

int mpoly_monomial_cmp_general(ulong * Aexp, flint_bitcnt_t Abits,
                   ulong * Bexp, flint_bitcnt_t Bbits, const mpoly_ctx_t mctx);

void mpoly_search_monomials(
                slong ** e_ind, ulong * e, slong * e_score,
                slong * t1, slong * t2, slong *t3,
                slong lower, slong upper,
                const ulong * a, slong a_len, const ulong * b, slong b_len,
                                               slong N, const ulong * cmpmask);

void mpoly_main_variable_split_LEX(slong * ind, ulong * pexp,
                                 const ulong * Aexp, slong l1, slong Alen,
                                  const ulong * mults, slong num, slong Abits);

void mpoly_main_variable_split_DEG(slong * ind, ulong * pexp,
                                  const ulong * Aexp,  slong l1, slong Alen,
                                            ulong deg, slong num, slong Abits);

int mpoly_term_exp_fits_si(ulong * exps, flint_bitcnt_t bits,
                                              slong n, const mpoly_ctx_t mctx);

int mpoly_term_exp_fits_ui(ulong * exps, flint_bitcnt_t bits,
                                              slong n, const mpoly_ctx_t mctx);

int mpoly_is_gen(ulong * exps, slong var,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);

int mpoly_monomials_valid_test(ulong * exps, slong len,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);

int mpoly_monomials_overflow_test(ulong * exps, slong len,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);

int mpoly_monomials_inorder_test(ulong * exps, slong len,
                                  flint_bitcnt_t bits, const mpoly_ctx_t mctx);

void mpoly_reverse(ulong * Aexp, const ulong * Bexp, slong len, slong N);

void mpoly_monomials_deflation(fmpz * shift, fmpz * stride,
                        const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                                                       const mpoly_ctx_t mctx);

void mpoly_monomials_deflate(ulong * Aexps, flint_bitcnt_t Abits,
                        const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength,
              const fmpz * shift, const fmpz * stride, const mpoly_ctx_t mctx);

void mpoly_monomials_inflate(ulong * Aexps, flint_bitcnt_t Abits,
                        const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength,
              const fmpz * shift, const fmpz * stride, const mpoly_ctx_t mctx);

void _mpoly_gen_shift_right(ulong * Aexp, flint_bitcnt_t Abits,
               slong Alength, slong var, ulong amount, const mpoly_ctx_t mctx);

void _mpoly_gen_shift_right_fmpz(ulong * Aexp, flint_bitcnt_t Abits,
        slong Alength, slong var, const fmpz_t amount, const mpoly_ctx_t mctx);

void _mpoly_gen_shift_left(ulong * Aexp, flint_bitcnt_t Abits,
               slong Alength, slong var, ulong amount, const mpoly_ctx_t mctx);

void mpoly_monomials_shift_right_ui(ulong * Aexps, flint_bitcnt_t Abits,
               slong Alength, const ulong * user_exps, const mpoly_ctx_t mctx);

void mpoly_monomials_shift_right_ffmpz(ulong * Aexps, flint_bitcnt_t Abits,
                slong Alength, const fmpz * user_exps, const mpoly_ctx_t mctx);

void mpoly1_fill_marks(ulong ** Dcoeffs, slong * Dlen, slong * Dalloc,
                        const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
                                                       const mpoly_ctx_t mctx);

void mpoly2_fill_marks(ulong ** Dcoeffs, slong * Dlen, slong * Dalloc,
                        const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
                                                       const mpoly_ctx_t mctx);

void mpoly_to_mpolyl_perm_deflate(
    ulong * Aexps,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t Actx,
    ulong * Bexps,
    flint_bitcnt_t Bbits,
    const mpoly_ctx_t Bctx,
    slong length,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

void mpoly_from_mpolyl_perm_inflate(
    ulong * Bexps,
    flint_bitcnt_t Bbits,
    const mpoly_ctx_t Bctx,
    ulong * Aexps,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t Actx,
    slong length,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

/* gcd ***********************************************************************/

#define MPOLY_GCD_USE_HENSEL  1
#define MPOLY_GCD_USE_BROWN   2
#define MPOLY_GCD_USE_ZIPPEL  4
#define MPOLY_GCD_USE_ZIPPEL2 8
#define MPOLY_GCD_USE_PRS     16
#define MPOLY_GCD_USE_ALL     31

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
    slong Adeflate_tdeg;
    slong Bdeflate_tdeg;

    double Adensity;
    double Bdensity;

    double hensel_time, brown_time, zippel_time, zippel2_time;
    slong * hensel_perm, * brown_perm, * zippel_perm, * zippel2_perm;
    unsigned int can_use;
    int Gdeflate_deg_bounds_are_nice; /* all of Gdeflate_deg_bound came from real gcd computations */

    char * data;
} mpoly_gcd_info_struct;

typedef mpoly_gcd_info_struct mpoly_gcd_info_t[1];

void mpoly_gcd_info_init(mpoly_gcd_info_t I, slong nvars);

void mpoly_gcd_info_clear(mpoly_gcd_info_t I);

void mpoly_gcd_info_limits(ulong * Amax_exp, ulong * Amin_exp,
                       slong * Amax_exp_count, slong * Amin_exp_count,
                       const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                                                       const mpoly_ctx_t mctx);
void mpoly_gcd_info_stride(ulong * strides,
          const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                             const ulong * Amax_exp, const ulong * Amin_exp,
          const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength,
                             const ulong * Bmax_exp, const ulong * Bmin_exp,
                                                       const mpoly_ctx_t mctx);

void mpoly_gcd_info_set_perm(mpoly_gcd_info_t I,
                         slong Alength, slong Blength, const mpoly_ctx_t mctx);

slong mpoly_gcd_info_get_brown_upper_limit(const mpoly_gcd_info_t I,
                                                       slong var, slong bound);

void mpoly_gcd_info_measure_hensel(mpoly_gcd_info_t I,
                         slong Alength, slong Blength, const mpoly_ctx_t mctx);

void mpoly_gcd_info_measure_brown(mpoly_gcd_info_t I,
                         slong Alength, slong Blength, const mpoly_ctx_t mctx);

void mpoly_gcd_info_measure_bma(mpoly_gcd_info_t I,
                         slong Alength, slong Blength, const mpoly_ctx_t mctx);

void mpoly_gcd_info_measure_zippel(mpoly_gcd_info_t I,
                         slong Alength, slong Blength, const mpoly_ctx_t mctx);

void mpoly_gcd_info_measure_zippel2(mpoly_gcd_info_t I,
                         slong Alength, slong Blength, const mpoly_ctx_t mctx);

int mpoly_monomial_cofactors(fmpz * Abarexps, fmpz * Bbarexps,
                                    const ulong * Aexps, flint_bitcnt_t Abits,
                                    const ulong * Bexps, flint_bitcnt_t Bbits,
                                        slong length,  const mpoly_ctx_t mctx);

/* factoring ****************************************************************/

#define MPOLY_FACTOR_USE_ZAS  1
#define MPOLY_FACTOR_USE_WANG 2
#define MPOLY_FACTOR_USE_ZIP  4
#define MPOLY_FACTOR_USE_ALL  7

int mpoly_is_proved_not_square(const ulong * Aexps,
                         slong Alen, flint_bitcnt_t Abits, slong N, ulong * t);

void mpoly_remove_var_powers(fmpz * var_powers, ulong * Aexps,
                     flint_bitcnt_t Abits, slong Alen, const mpoly_ctx_t mctx);

slong _mpoly_compress_exps(slong * V, slong * D, slong * deg,
                                                  slong * S, slong n, slong l);

int mpoly_test_irreducible(ulong * Aexps, flint_bitcnt_t Abits,
                                            slong Alen, const mpoly_ctx_t ctx);

int _mpoly_test_irreducible(slong * Aexps, slong stride, slong Alen,
                            slong nvars, flint_rand_t state, slong tries_left);

typedef struct {
    slong mvars;
    slong nvars;
    slong * exps;
    slong exps_alloc;
    slong * rest;
    slong rest_alloc;
    slong * umat;
    slong * deltas;
    slong * degs;
    int is_trivial;
    int is_perm;
    int is_irred;
} mpoly_compression_struct;

typedef mpoly_compression_struct mpoly_compression_t[1];

void mpoly_compression_init(mpoly_compression_t M);

void mpoly_compression_clear(mpoly_compression_t M);

void mpoly_compression_set(mpoly_compression_t M, const ulong * Aexps,
                     flint_bitcnt_t Abits, slong Alen, const mpoly_ctx_t mctx);

void mpoly_bivar_cld_bounds(slong * l, slong n);

MPOLY_INLINE
void _slong_array_fit_length(slong ** array, slong * alloc, slong len)
{
    if (len <= *alloc)
        return;
    len = FLINT_MAX(len, *alloc + *alloc/4 + 1);
    *alloc = len;
    *array = (slong *) flint_realloc(*array, len*sizeof(slong));
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

void * _mpoly_heap_pop1(mpoly_heap1_s * heap, slong * heap_len, ulong maskhi);

void _mpoly_heap_insert1(mpoly_heap1_s * heap, ulong exp, void * x,
                              slong * next_loc, slong * heap_len, ulong maskhi);

void * _mpoly_heap_pop(mpoly_heap_s * heap, slong * heap_len, slong N,
                                                         const ulong * cmpmask);

int _mpoly_heap_insert(mpoly_heap_s * heap, ulong * exp, void * x,
       slong * next_loc, slong * heap_len, slong N, const ulong * cmpmask);

/* generic parsing ***********************************************************/

typedef struct {
    char * coeffs;
    fmpz * exps;
    slong length;
    slong alloc;
} mpoly_univar_struct;

typedef mpoly_univar_struct mpoly_univar_t[1];

typedef struct {
    slong elem_size;
    const void * ctx;
    void (*init)(void *, const void *);
    void (*clear)(void *, const void *);
    int (*is_zero)(const void *, const void *);
    void (*zero)(void *, const void *);
    void (*one)(void *, const void *);
    void (*set_fmpz)(void *, const fmpz_t, const void *);
    void (*set)(void *, const void *, const void *);
    void (*swap)(void *, void *, const void *);
    void (*neg)(void *, const void *, const void *);
    void (*add)(void *, const void *, const void *, const void *);
    void (*sub)(void *, const void *, const void *, const void *);
    void (*mul_fmpz)(void *, const void *, const fmpz_t, const void *);
    void (*mul)(void *, const void *, const void *, const void *);
    void (*divexact)(void *, const void *, const void *, const void *);
    int (*divides)(void *, const void *, const void *, const void *);
    int (*pow_fmpz)(void *, const void *, const fmpz_t, const void *);
    slong (*length)(const void *, const void *);
} mpoly_void_ring_t[1];

void * mpoly_void_ring_elem_init(mpoly_void_ring_t R);

void mpoly_void_ring_elem_clear(void * a, mpoly_void_ring_t R);

void mpoly_univar_init(mpoly_univar_t A, mpoly_void_ring_t R);

void mpoly_univar_clear(mpoly_univar_t A, mpoly_void_ring_t R);

void mpoly_univar_swap(mpoly_univar_t A, mpoly_univar_t B);

void mpoly_univar_fit_length(mpoly_univar_t A, slong len,
                                                          mpoly_void_ring_t R);

void mpoly_univar_init2(mpoly_univar_t A, slong len,
                                                          mpoly_void_ring_t R);

int mpoly_univar_pseudo_gcd_ducos(mpoly_univar_t G,
                      mpoly_univar_t B, mpoly_univar_t A, mpoly_void_ring_t R);

int mpoly_univar_resultant(void * r, mpoly_univar_t fx,
                                       mpoly_univar_t gx, mpoly_void_ring_t R);

int mpoly_univar_discriminant(void * d, mpoly_univar_t fx,
                                                          mpoly_void_ring_t R);

typedef struct {
    char * str;
    slong str_len;
} string_with_length_struct;

typedef struct {
    mpoly_void_ring_t R;
    slong * stack;
    slong stack_len;
    slong stack_alloc;
    char * estore;
    slong estore_len;
    slong estore_alloc;
    void * tmp;
    string_with_length_struct * terminal_strings;
    char * terminal_values;
    slong terminals_alloc;
    slong terminals_len;
} mpoly_parse_struct;

typedef mpoly_parse_struct mpoly_parse_t[1];

void mpoly_parse_init(mpoly_parse_t E);

void mpoly_parse_clear(mpoly_parse_t E);

void mpoly_parse_add_terminal(mpoly_parse_t E,
                                               const char * s, const void * v);

int mpoly_parse_parse(mpoly_parse_t E, void * res,
                                                    const char * s, slong len);

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

