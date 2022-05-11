/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_H
#define NMOD_H

#ifdef NMOD_INLINES_C
#define NMOD_INLINE FLINT_DLL
#else
#define NMOD_INLINE static __inline__
#endif

#include "ulong_extras_mini.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define NMOD_RED2(r, a_hi, a_lo, mod)                                   \
    do                                                                  \
    {                                                                   \
        ulong q0xx, q1xx, r1xx;                                     \
        const ulong u1xx = ((a_hi)<<(mod).norm)                     \
            + ((mod.norm == 0)                                          \
                    ? UWORD(0)                                          \
                    : (ulong) (a_lo) >> (FLINT_BITS - (mod).norm)); \
        const ulong u0xx = ((a_lo)<<(mod).norm);                    \
        const ulong nxx = ((mod).n<<(mod).norm);                    \
        umul_ppmm(q1xx, q0xx, (mod).ninv, u1xx);                        \
        add_ssaaaa(q1xx, q0xx, q1xx, q0xx, u1xx, u0xx);                 \
        r1xx = (u0xx - (q1xx + 1)*nxx);                                 \
        if (r1xx > q0xx)                                                \
            r1xx += nxx;                                                \
        if (r1xx < nxx)                                                 \
            r = (r1xx>>(mod).norm);                                     \
        else                                                            \
            r = ((r1xx - nxx)>>(mod).norm);                             \
    } while (0)

#define NMOD_RED(r, a, mod)             \
   do {                                 \
      NMOD_RED2(r, UWORD(0), a, mod);   \
   } while (0)

#define NMOD2_RED2(r, a_hi, a_lo, mod)  \
    do {                                \
       ulong v_hi;	                \
       NMOD_RED(v_hi, a_hi, mod);       \
       NMOD_RED2(r, v_hi, a_lo, mod);   \
    } while (0)

#define NMOD_RED3(r, a_hi, a_me, a_lo, mod) \
    do {                                    \
       ulong v_hi;	                    \
       NMOD_RED2(v_hi, a_hi, a_me, mod);    \
       NMOD_RED2(r, v_hi, a_lo, mod);       \
    } while (0)

#define NMOD_BITS(mod) (FLINT_BITS - ((mod).norm))
#define NMOD_CAN_USE_SHOUP(mod) ((mod).norm > 0)

#define NMOD_MUL_PRENORM(res, a, b, mod) \
    do { \
        ulong q0xx, q1xx, rxx, p_hixx, p_loxx; \
        ulong nxx, ninvxx; \
        unsigned int normxx; \
        ninvxx = (mod).ninv; \
        normxx = (mod).norm; \
        nxx = (mod).n << normxx; \
        umul_ppmm(p_hixx, p_loxx, (a), (b)); \
        umul_ppmm(q1xx, q0xx, ninvxx, p_hixx); \
        add_ssaaaa(q1xx, q0xx, q1xx, q0xx, p_hixx, p_loxx); \
        rxx = (p_loxx - (q1xx + 1) * nxx); \
        if (rxx > q0xx) \
            rxx += nxx; \
        rxx = (rxx < nxx ? rxx : rxx - nxx) >> normxx; \
        (res) = rxx; \
    } while (0)

#define NMOD_MUL_FULLWORD(res, a, b, mod) \
    do { \
        ulong q0xx, q1xx, rxx, p_hixx, p_loxx; \
        ulong nxx, ninvxx; \
        ninvxx = (mod).ninv; \
        nxx = (mod).n; \
        umul_ppmm(p_hixx, p_loxx, (a), (b)); \
        umul_ppmm(q1xx, q0xx, ninvxx, p_hixx); \
        add_ssaaaa(q1xx, q0xx, q1xx, q0xx, p_hixx, p_loxx); \
        rxx = (p_loxx - (q1xx + 1) * nxx); \
        if (rxx > q0xx) \
            rxx += nxx; \
        rxx = (rxx < nxx ? rxx : rxx - nxx); \
        (res) = rxx; \
    } while (0)

NMOD_INLINE
ulong _nmod_add(ulong a, ulong b, nmod_t mod)
{
   const ulong sum = a + b;
   return sum - mod.n + ((((slong)(sum - mod.n))>>(FLINT_BITS - 1)) & mod.n);
}

NMOD_INLINE
ulong _nmod_sub(ulong a, ulong b, nmod_t mod)
{
   const ulong diff = a - b;
   return  ((((slong)diff)>>(FLINT_BITS - 1)) & mod.n) + diff;
}

NMOD_INLINE
ulong nmod_add(ulong a, ulong b, nmod_t mod)
{
   const ulong neg = mod.n - a;
   if (neg > b)
      return a + b;
   else 
      return b - neg;
}

NMOD_INLINE
ulong nmod_sub(ulong a, ulong b, nmod_t mod)
{
   const ulong diff = a - b;
   
   if (a < b)
      return mod.n + diff;
   else
      return diff;
}

NMOD_INLINE
ulong nmod_neg(ulong a, nmod_t mod)
{
   if (a)
      return mod.n - a;
   else
      return 0;
}

NMOD_INLINE
ulong nmod_mul(ulong a, ulong b, nmod_t mod)
{
    ulong res;
    NMOD_MUL_PRENORM(res, a, b << mod.norm, mod);
    return res;
}

NMOD_INLINE
ulong _nmod_mul_fullword(ulong a, ulong b, nmod_t mod)
{
    ulong res;
    NMOD_MUL_FULLWORD(res, a, b, mod);
    return res;
}

NMOD_INLINE
ulong nmod_addmul(ulong a, ulong b, ulong c, nmod_t mod)
{
    return nmod_add(a, nmod_mul(b, c, mod), mod);
}

#define NMOD_ADDMUL(r, a, b, mod) \
    do { \
       (r) = nmod_addmul((r), (a), (b), (mod)); \
    } while (0)

NMOD_INLINE
ulong nmod_inv(ulong a, nmod_t mod)
{
    return n_invmod(a, mod.n);
}

NMOD_INLINE
ulong nmod_div(ulong a, ulong b, nmod_t mod)
{
    return nmod_mul(a, n_invmod(b, mod.n), mod);
}

NMOD_INLINE
ulong nmod_pow_ui(ulong a, ulong exp, nmod_t mod)
{
    return n_powmod2_ui_preinv(a, exp, mod.n, mod.ninv);
}

NMOD_INLINE ulong nmod_pow_fmpz(ulong a, const fmpz_t exp, nmod_t mod)
{
    return n_powmod2_fmpz_preinv(a, exp, mod.n, mod.ninv);
}


NMOD_INLINE
void nmod_init(nmod_t * mod, ulong n)
{
   mod->n = n;
   mod->ninv = n_preinvert_limb(n);
   count_leading_zeros(mod->norm, n);
}

#ifdef __cplusplus
}
#endif

#endif
