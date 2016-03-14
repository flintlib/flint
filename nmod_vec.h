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

    Copyright (C) 2010 William Hart

******************************************************************************/

#ifndef NMOD_VEC_H
#define NMOD_VEC_H

#ifdef NMOD_VEC_INLINES_C
#define NMOD_VEC_INLINE FLINT_DLL
#else
#define NMOD_VEC_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "longlong.h"
#include "ulong_extras.h"
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
   mp_limb_t n;
   mp_limb_t ninv;
   mp_bitcnt_t norm;
} nmod_t;


#define NMOD_VEC_NORM(vec, i)                   \
do {                                            \
    while ((i) && vec[(i) - 1] == UWORD(0))     \
        (i)--;                                  \
} while (0)

#define NMOD_RED2(r, a_hi, a_lo, mod) \
   do { \
      mp_limb_t q0, q1, r1; \
      const mp_limb_t u1 = ((a_hi)<<(mod).norm) + r_shift((a_lo), FLINT_BITS - (mod).norm);	\
      const mp_limb_t u0 = ((a_lo)<<(mod).norm); \
      const mp_limb_t nxx = ((mod).n<<(mod).norm); \
      umul_ppmm(q1, q0, (mod).ninv, u1); \
      add_ssaaaa(q1, q0, q1, q0, u1, u0); \
      r1 = (u0 - (q1 + 1)*nxx); \
      if (r1 > q0) r1 += nxx; \
      if (r1 < nxx) r = (r1>>(mod).norm); \
      else r = ((r1 - nxx)>>(mod).norm); \
   } while (0)

#define NMOD_RED(r, a, mod) \
   do { \
      NMOD_RED2(r, 0, a, mod); \
   } while (0)

#define NMOD2_RED2(r, a_hi, a_lo, mod) \
    do { \
       mp_limb_t v_hi;	\
       NMOD_RED(v_hi, a_hi, mod); \
       NMOD_RED2(r, v_hi, a_lo, mod); \
    } while (0)

#define NMOD_RED3(r, a_hi, a_me, a_lo, mod) \
    do { \
       mp_limb_t v_hi;	\
       NMOD_RED2(v_hi, a_hi, a_me, mod); \
       NMOD_RED2(r, v_hi, a_lo, mod); \
    } while (0)

#define NMOD_ADDMUL(r, a, b, mod) \
    do { \
       mp_limb_t a_hi, a_lo; \
       umul_ppmm(a_hi, a_lo, a, b); \
       add_ssaaaa(a_hi, a_lo, a_hi, a_lo, (mp_limb_t) 0, r); \
       NMOD_RED2(r, a_hi, a_lo, mod); \
    } while (0)

NMOD_VEC_INLINE
mp_limb_t _nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
   const mp_limb_t sum = a + b;
   return sum - mod.n + ((((mp_limb_signed_t)(sum - mod.n))>>(FLINT_BITS - 1)) & mod.n);
}

NMOD_VEC_INLINE
mp_limb_t _nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
   const mp_limb_t diff = a - b;
   return  ((((mp_limb_signed_t)diff)>>(FLINT_BITS - 1)) & mod.n) + diff;
}

NMOD_VEC_INLINE
mp_limb_t nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
   const mp_limb_t neg = mod.n - a;
   if (neg > b)
      return a + b;
   else 
      return b - neg;
}

NMOD_VEC_INLINE
mp_limb_t nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
   const mp_limb_t diff = a - b;
   
   if (a < b)
      return mod.n + diff;
   else
      return diff;
}

NMOD_VEC_INLINE
mp_limb_t nmod_neg(mp_limb_t a, nmod_t mod)
{
   if (a)
      return mod.n - a;
   else
      return 0;
}

NMOD_VEC_INLINE
mp_limb_t nmod_mul(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
    return n_mulmod2_preinv(a, b, mod.n, mod.ninv);
}

NMOD_VEC_INLINE
mp_limb_t nmod_inv(mp_limb_t a, nmod_t mod)
{
    return n_invmod(a, mod.n);
}

NMOD_VEC_INLINE
mp_limb_t nmod_div(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
    b = n_invmod(b, mod.n);
    return n_mulmod2_preinv(a, b, mod.n, mod.ninv);
}

NMOD_VEC_INLINE
mp_limb_t nmod_pow_ui(mp_limb_t a, ulong exp, nmod_t mod)
{
    return n_powmod2_ui_preinv(a, exp, mod.n, mod.ninv);
}

NMOD_VEC_INLINE
void nmod_init(nmod_t * mod, mp_limb_t n)
{
   mod->n = n;
   mod->ninv = n_preinvert_limb(n);
   count_leading_zeros(mod->norm, n);
}

NMOD_VEC_INLINE
mp_ptr _nmod_vec_init(slong len)
{
   return (mp_ptr) flint_malloc(len * sizeof(mp_limb_t));
}

NMOD_VEC_INLINE
void _nmod_vec_clear(mp_ptr vec)
{
   flint_free(vec);
}

FLINT_DLL void _nmod_vec_randtest(mp_ptr vec, flint_rand_t state, slong len, nmod_t mod);

NMOD_VEC_INLINE
void _nmod_vec_zero(mp_ptr vec, slong len)
{
   flint_mpn_zero(vec, len);
}

FLINT_DLL mp_bitcnt_t _nmod_vec_max_bits(mp_srcptr vec, slong len);

NMOD_VEC_INLINE
void _nmod_vec_set(mp_ptr res, mp_srcptr vec, slong len)
{
   flint_mpn_copyi(res, vec, len);
}

NMOD_VEC_INLINE
void _nmod_vec_swap(mp_ptr a, mp_ptr b, slong length)
{
    slong i;
    for (i = 0; i < length; i++)
    {
        mp_limb_t t = a[i];
        a[i] = b[i];
        b[i] = t;
    }
}

NMOD_VEC_INLINE
int _nmod_vec_equal(mp_srcptr vec, mp_srcptr vec2, slong len)
{
   slong i;

   for (i = 0; i < len; i++)
      if (vec[i] != vec2[i]) return 0;

   return 1;
}

NMOD_VEC_INLINE
int _nmod_vec_is_zero(mp_srcptr vec, slong len)
{
   slong i;

   for (i = 0; i < len; i++)
      if (vec[i] != 0) return 0;

   return 1;
}

FLINT_DLL void _nmod_vec_reduce(mp_ptr res, mp_srcptr vec, 
                                        slong len, nmod_t mod);

FLINT_DLL void _nmod_vec_add(mp_ptr res, mp_srcptr vec1, 
                        mp_srcptr vec2, slong len, nmod_t mod);

FLINT_DLL void _nmod_vec_sub(mp_ptr res, mp_srcptr vec1, 
                        mp_srcptr vec2, slong len, nmod_t mod);

FLINT_DLL void _nmod_vec_neg(mp_ptr res, mp_srcptr vec, 
                                            slong len, nmod_t mod);

FLINT_DLL void _nmod_vec_scalar_mul_nmod(mp_ptr res, mp_srcptr vec, 
                            slong len, mp_limb_t c, nmod_t mod);

FLINT_DLL void _nmod_vec_scalar_mul_nmod_shoup(mp_ptr res, mp_srcptr vec, 
                            slong len, mp_limb_t c, nmod_t mod);

FLINT_DLL void _nmod_vec_scalar_addmul_nmod(mp_ptr res, mp_srcptr vec, 
                            slong len, mp_limb_t c, nmod_t mod);

FLINT_DLL int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod);


#define NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, nlimbs)                \
    do                                                                      \
    {                                                                       \
        mp_limb_t s0, s1, s2, t0, t1;                                       \
        s0 = s1 = s2 = UWORD(0);                                                 \
        switch (nlimbs)                                                     \
        {                                                                   \
            case 1:                                                         \
                for (i = 0; i < len; i++)                                   \
                {                                                           \
                    s0 += (expr1) * (expr2);                                \
                }                                                           \
                NMOD_RED(s0, s0, mod);                                      \
                break;                                                      \
            case 2:                                                         \
                if (mod.n <= (UWORD(1) << (FLINT_BITS / 2)))                     \
                {                                                           \
                    for (i = 0; i < len; i++)                               \
                    {                                                       \
                        t0 = (expr1) * (expr2);                             \
                        add_ssaaaa(s1, s0, s1, s0, 0, t0);                  \
                    }                                                       \
                }                                                           \
                else                                                        \
                {                                                           \
                    for (i = 0; i < len; i++)                               \
                    {                                                       \
                        umul_ppmm(t1, t0, (expr1), (expr2));                \
                        add_ssaaaa(s1, s0, s1, s0, t1, t0);                 \
                    }                                                       \
                }                                                           \
                NMOD2_RED2(s0, s1, s0, mod);                                \
                break;                                                      \
            default:                                                        \
                for (i = 0; i < len; i++)                                   \
                {                                                           \
                    umul_ppmm(t1, t0, (expr1), (expr2));                    \
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);       \
                }                                                           \
                NMOD_RED(s2, s2, mod);                                      \
                NMOD_RED3(s0, s2, s1, s0, mod);                             \
                break;                                                      \
        }                                                                   \
        res = s0;                                                           \
    } while (0);

FLINT_DLL mp_limb_t _nmod_vec_dot(mp_srcptr vec1, mp_srcptr vec2,
    slong len, nmod_t mod, int nlimbs);

FLINT_DLL mp_limb_t _nmod_vec_dot_ptr(mp_srcptr vec1, const mp_ptr * vec2, slong offset,
    slong len, nmod_t mod, int nlimbs);

#ifdef __cplusplus
}
#endif

#endif

