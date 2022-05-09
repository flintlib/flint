/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPN_EXTRAS_H
#define MPN_EXTRAS_H

#ifdef MPN_EXTRAS_INLINES_C
#define MPN_EXTRAS_INLINE FLINT_DLL
#else
#define MPN_EXTRAS_INLINE static __inline__
#endif

#include "gmp.h"
#include "flint.h"
#ifdef LONGSLONG
# define flint_mpz_divisible_ui_p mpz_divisible_ui_p
#else
# include "gmpcompat.h"
#endif

#ifdef __cplusplus
 extern "C" {
#endif

#define MPN_NORM(a, an)                         \
    do {                                        \
        while ((an) != 0 && (a)[(an) - 1] == 0) \
           (an)--;                              \
    } while (0)

#define MPN_SWAP(a, an, b, bn) \
    do {                       \
        ulong_ptr __t;            \
        mp_size_t __tn;        \
        __t = (a);             \
        (a) = (b);             \
        (b) = __t;             \
        __tn = (an);           \
        (an) = (bn);           \
        (bn) = __tn;           \
    } while (0)

/*
    return the high limb of a two limb left shift by n < GMP_LIMB_BITS bits.
    Note: if GMP_NAIL_BITS != 0, the rest of flint is already broken anyways.
*/
#define MPN_LEFT_SHIFT_HI(hi, lo, n)                                \
    ((n) > 0 ? (((hi) << (n)) | ((lo) >> (GMP_LIMB_BITS - (n))))    \
             : (hi))

#define MPN_RIGHT_SHIFT_LOW(hi, lo, n)                                \
    ((n) > 0 ? (((lo) >> (n)) | ((hi) << (GMP_LIMB_BITS - (n))))    \
             : (lo))


/* Not defined in gmp.h
ulong  __gmpn_modexact_1_odd(ulong_srcptr src, mp_size_t size,
                                 ulong divisor);
#define mpn_modexact_1_odd __gmpn_modexact_1_odd
*/

#ifdef mpn_modexact_1_odd
#define flint_mpn_divisible_1_p(x, xsize, d) (mpn_modexact_1_odd(x, xsize, d) == 0)
#else
MPN_EXTRAS_INLINE int
flint_mpn_divisible_1_p(ulong_srcptr x, mp_size_t xsize, ulong d)
{
    __mpz_struct s;
    s._mp_size = xsize;
    s._mp_d = (ulong_ptr) x;
    return flint_mpz_divisible_ui_p(&s, d);
}
#endif

MPN_EXTRAS_INLINE
int flint_mpn_zero_p(ulong_srcptr x, mp_size_t xsize)
{
    slong i;
    for (i = 0; i < xsize; i++)
    {
        if (x[i])
            return 0;
    }
    return 1;
}

MPN_EXTRAS_INLINE
mp_size_t flint_mpn_divexact_1(ulong_ptr x, mp_size_t xsize, ulong d)
{
    mpn_divrem_1(x, 0, x, xsize, d);
    if (x[xsize - 1] == UWORD(0))
        xsize -= 1;
    return xsize;
}

#if defined(__MPIR_VERSION)

#if !defined(__MPIR_RELEASE ) || __MPIR_RELEASE < 20600
#define mpn_sumdiff_n __MPN(sumdiff_n)
extern
ulong mpn_sumdiff_n(ulong_ptr, ulong_ptr, ulong_srcptr, ulong_srcptr, mp_mock_size_t);
#endif

#else

/* TODO: Don't inline this. Else, we need to include flint-impl.h, which the
 * user should not see. */
#include "flint-impl.h"
MPN_EXTRAS_INLINE
ulong mpn_sumdiff_n(ulong_ptr s, ulong_ptr d, ulong_srcptr x, ulong_srcptr y, mp_mock_size_t n)
{
    ulong ret;
    ulong_ptr t;

    if (n == 0)
        return 0;

    if ((s == x && d == y) || (s == y && d == x))
    {
        t = (ulong_ptr) flint_malloc(n * sizeof(ulong));
        ret = mpn_sub_n(t, x, y, n);
        ret += 2 * mpn_add_n(s, x, y, n);
        FLINT_MPN_COPYI(d, t, n);
        flint_free(t);
        return ret;
    }

    if (s == x || s == y)
    {
        ret = mpn_sub_n(d, x, y, n);
        ret += 2 * mpn_add_n(s, x, y, n);
        return ret;
    }

    ret = 2 * mpn_add_n(s, x, y, n);
    ret += mpn_sub_n(d, x, y, n);
    return ret;
}

#endif

/* I do not think this should be inlined as mpn_add_1 and mpn_sub_1 is usually
 * inlined. */
MPN_EXTRAS_INLINE
void mpn_addmod_2expp1_1(ulong * r, mp_mock_size_t limbs, slong c)
{
   ulong sum = r[0] + c;

   /* check if adding c would cause a carry to propagate */
   if ((slong)(sum ^ r[0]) >= 0)
      r[0] = sum;
   else
   {
      if (c >= 0)
          mpn_add_1(r, r, limbs + 1, c);
      else
          mpn_sub_1(r, r, limbs + 1, -c);
   }
}

FLINT_DLL void flint_mpn_debug(ulong_srcptr x, mp_size_t xsize);

FLINT_DLL mp_size_t flint_mpn_remove_2exp(ulong_ptr x, mp_size_t xsize,
		                                      flint_bitcnt_t *bits);

FLINT_DLL mp_size_t flint_mpn_remove_power_ascending(ulong_ptr x,
		    mp_size_t xsize, ulong_ptr p, mp_size_t psize, ulong *exp);

FLINT_DLL int flint_mpn_factor_trial(ulong_srcptr x, mp_size_t xsize,
		                                   slong start, slong stop);

FLINT_DLL int flint_mpn_factor_trial_tree(slong * factors,
                            ulong_srcptr x, mp_size_t xsize, slong num_primes);

FLINT_DLL mp_size_t flint_mpn_fmms1(ulong_ptr y, ulong a1, ulong_srcptr x1,
                                      ulong a2, ulong_srcptr x2, mp_size_t n);

FLINT_DLL int flint_mpn_divides(ulong_ptr q, ulong_srcptr array1, 
         mp_size_t limbs1, ulong_srcptr arrayg, mp_size_t limbsg, ulong_ptr temp);

FLINT_DLL mp_size_t flint_mpn_gcd_full2(ulong_ptr arrayg,
		                 ulong_srcptr array1, mp_size_t limbs1,
			   ulong_srcptr array2, mp_size_t limbs2, ulong_ptr temp);

FLINT_DLL mp_size_t flint_mpn_gcd_full(ulong_ptr arrayg, 
    ulong_srcptr array1, mp_size_t limbs1, ulong_srcptr array2, mp_size_t limbs2);

FLINT_DLL ulong flint_mpn_preinv1(ulong d, ulong d2);

FLINT_DLL ulong flint_mpn_divrem_preinv1(ulong_ptr q, ulong_ptr a, 
           mp_size_t m, ulong_srcptr b, mp_size_t n, ulong dinv);

#define flint_mpn_divrem21_preinv(q, a_hi, a_lo, dinv) \
   do { \
      ulong __q2, __q3, __q4; \
      umul_ppmm((q), __q2, (a_hi), (dinv)); \
      umul_ppmm(__q3, __q4, (a_lo), (dinv)); \
      add_ssaaaa((q), __q2, (q), __q2, 0, __q3); \
      add_ssaaaa((q), __q2, (q), __q2, (a_hi), (a_lo)); \
   } while (0)

FLINT_DLL void flint_mpn_mulmod_preinv1(ulong_ptr r, 
        ulong_srcptr a, ulong_srcptr b, mp_size_t n, 
        ulong_srcptr d, ulong dinv, ulong norm);

FLINT_DLL void flint_mpn_preinvn(ulong_ptr dinv, ulong_srcptr d, mp_size_t n);

FLINT_DLL void flint_mpn_mod_preinvn(ulong_ptr r, ulong_srcptr a, mp_size_t m, 
                                     ulong_srcptr d, mp_size_t n, ulong_srcptr dinv);

FLINT_DLL ulong flint_mpn_divrem_preinvn(ulong_ptr q, ulong_ptr r, ulong_srcptr a, mp_size_t m, 
                                     ulong_srcptr d, mp_size_t n, ulong_srcptr dinv);

FLINT_DLL void flint_mpn_mulmod_preinvn(ulong_ptr r, 
        ulong_srcptr a, ulong_srcptr b, mp_size_t n, 
        ulong_srcptr d, ulong_srcptr dinv, ulong norm);

FLINT_DLL int flint_mpn_mulmod_2expp1_basecase(ulong_ptr xp, ulong_srcptr yp, ulong_srcptr zp, 
    int c, flint_bitcnt_t b, ulong_ptr tp);

MPN_EXTRAS_INLINE
void flint_mpn_rrandom(ulong *rp, gmp_randstate_t state, mp_size_t n)
{
  __mpz_struct str;
  str._mp_d = rp;
  str._mp_alloc = n;
  str._mp_size =n;
  mpz_rrandomb(&str,state,FLINT_BITS*n);
}

MPN_EXTRAS_INLINE
void flint_mpn_urandomb(ulong *rp, gmp_randstate_t state, flint_bitcnt_t n)
{
  __mpz_struct str;
  str._mp_d = rp;
  str._mp_alloc = (n + FLINT_BITS - 1)/FLINT_BITS;
  str._mp_size = (n + FLINT_BITS - 1)/FLINT_BITS;
  mpz_rrandomb(&str,state,n);
}


/******************************************************************************
    Divisions where the quotient is expected to be small. All function do:
        input: n > d > 0
        output: q = n/d, r = n%d
    for various small sizes of n and d.
    Not in a function because compiler refuses to inline eudiv_qrrnndd.
    Each macro takes a prefix t for its local vars.
******************************************************************************/

#define eudiv_qrnd(q, r, n, d, t)           \
do {                                        \
    ulong t##q, t##a = n, t##b = d;     \
                                            \
    FLINT_ASSERT(t##a > t##b);              \
    FLINT_ASSERT(t##b > 0);                 \
                                            \
    t##a -= t##b;                           \
    for (t##q = 1; t##q < 5; t##q++)        \
    {                                       \
        if (t##a < t##b)                    \
            goto t##quotient_found;         \
        t##a -= t##b;                       \
    }                                       \
    t##q += t##a / t##b;                    \
    t##a = t##a % t##b;                     \
                                            \
t##quotient_found:                          \
                                            \
    q = t##q;                               \
    r = t##a;                               \
} while (0)


#define eudiv_qqrnnd(q1, q0, r0, n1, n0, d0, t)         \
do {                                                    \
    ulong t##a1 = n1, t##a0 = n0, t##b0 = d0;       \
    ulong t##q1, t##q0, t##r0, t##u;                \
                                                        \
    FLINT_ASSERT(t##a1 > 0 || t##a0 > t##b0);           \
                                                        \
    udiv_qrnnd(t##q1, t##u, 0, t##a1, t##b0);           \
    udiv_qrnnd(t##q0, t##r0, t##u, t##a0, t##b0);       \
                                                        \
    q1 = t##q1;                                         \
    q0 = t##q0;                                         \
    r0 = t##r0;                                         \
} while (0)

/* d must be normalized, i.e. d1 != 0 */
#define eudiv_qrrnndd(q0, r1, r0, n1, n0, d1, d0, t)                        \
do {                                                                        \
    int t##i;                                                               \
    ulong t##a1 = n1, t##a0 = n0, t##b1 = d1, t##b0 = d0;               \
    ulong t##r1, t##r0, t##u2, t##u1, t##u0, t##q, t##v1, t##v0;        \
                                                                            \
    FLINT_ASSERT(t##a1 != 0);                                               \
    FLINT_ASSERT(t##b1 != 0);                                               \
    FLINT_ASSERT(t##b1 < t##a1 || (t##b1 == t##a1 && t##b0 < t##a0));       \
                                                                            \
    t##q = 1;                                                               \
                                                                            \
    sub_ddmmss(t##r1,t##r0, t##a1,t##a0, t##b1,t##b0);                      \
                                                                            \
t##subtract:                                                                \
                                                                            \
    for (t##i = 2; t##i <= 4; t##i++)                                       \
    {                                                                       \
        sub_dddmmmsss(t##u2,t##u1,t##u0, 0,t##r1,t##r0, 0,t##b1,t##b0);     \
        if (t##u2 != 0)                                                     \
            goto t##quotient_found;                                         \
        t##q += 1;                                                          \
        t##r0 = t##u0;                                                      \
        t##r1 = t##u1;                                                      \
    }                                                                       \
                                                                            \
    if (t##r1 != 0)                                                         \
    {                                                                       \
        int t##ncnt, t##dcnt;                                               \
        ulong t##qq = 0;                                                \
                                                                            \
        count_leading_zeros(t##ncnt, t##r1);                                \
        count_leading_zeros(t##dcnt, t##b1);                                \
        t##dcnt -= t##ncnt;                                                 \
        if (t##dcnt <= 0)                                                   \
            goto t##subtract;                                               \
                                                                            \
        t##v1 = (t##b1 << t##dcnt) | (t##b0 >> (FLINT_BITS - t##dcnt));     \
        t##v0 = t##b0 << t##dcnt;                                           \
                                                                            \
        do {                                                                \
            sub_dddmmmsss(t##u2,t##u1,t##u0, 0,t##r1,t##r0, 0,t##v1,t##v0); \
            t##qq = 2*t##qq + 1 + t##u2;                                    \
            t##r1 = t##u2 ? t##r1 : t##u1;                                  \
            t##r0 = t##u2 ? t##r0 : t##u0;                                  \
            t##v0 = (t##v1 << (FLINT_BITS - 1)) | (t##v0 >> 1);             \
            t##v1 = t##v1 >> 1;                                             \
            t##dcnt--;                                                      \
        } while (t##dcnt >= 0);                                             \
                                                                            \
        t##q += t##qq;                                                      \
    }                                                                       \
                                                                            \
t##quotient_found:                                                          \
                                                                            \
    FLINT_ASSERT(t##r1 < t##b1 || (t##r1 == t##b1 && t##r0 < t##b0));       \
                                                                            \
    q0 = t##q;                                                              \
    r0 = t##r0;                                                             \
    r1 = t##r1;                                                             \
} while (0)


#ifdef __cplusplus
}
#endif

#endif
