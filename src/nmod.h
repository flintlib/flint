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
#define NMOD_INLINE
#else
#define NMOD_INLINE static inline
#endif

#include "ulong_extras.h"
#include "nmod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NMOD_RED2(r, a_hi, a_lo, mod) \
   do { \
      mp_limb_t q0xx, q1xx, r1xx; \
      const mp_limb_t u1xx = ((a_hi)<<(mod).norm) + r_shift((a_lo), FLINT_BITS - (mod).norm);	\
      const mp_limb_t u0xx = ((a_lo)<<(mod).norm); \
      const mp_limb_t nxx = ((mod).n<<(mod).norm); \
      umul_ppmm(q1xx, q0xx, (mod).ninv, u1xx); \
      add_ssaaaa(q1xx, q0xx, q1xx, q0xx, u1xx, u0xx); \
      r1xx = (u0xx - (q1xx + 1)*nxx); \
      if (r1xx > q0xx) r1xx += nxx; \
      if (r1xx < nxx) r = (r1xx>>(mod).norm); \
      else r = ((r1xx - nxx)>>(mod).norm); \
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

#define NMOD_BITS(mod) (FLINT_BITS - ((mod).norm))
#define NMOD_CAN_USE_SHOUP(mod) ((mod).norm > 0)

#define NMOD_MUL_PRENORM(res, a, b, mod) \
    do { \
        mp_limb_t q0xx, q1xx, rxx, p_hixx, p_loxx; \
        mp_limb_t nxx, ninvxx; \
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
        mp_limb_t q0xx, q1xx, rxx, p_hixx, p_loxx; \
        mp_limb_t nxx, ninvxx; \
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

NMOD_INLINE mp_limb_t nmod_set_ui(ulong x, nmod_t mod)
{
    if (x < mod.n)
        return x;

    NMOD_RED(x, x, mod);
    return x;
}

NMOD_INLINE
mp_limb_t nmod_set_si(slong x, nmod_t mod)
{
    ulong res = FLINT_ABS(x);
    NMOD_RED(res, res, mod);
    return (res == 0 || x > 0) ? res : mod.n - res;
}

NMOD_INLINE
mp_limb_t _nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
   const mp_limb_t sum = a + b;
   return sum - mod.n + ((((mp_limb_signed_t)(sum - mod.n))>>(FLINT_BITS - 1)) & mod.n);
}

NMOD_INLINE
mp_limb_t _nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
   const mp_limb_t diff = a - b;
   return  ((((mp_limb_signed_t)diff)>>(FLINT_BITS - 1)) & mod.n) + diff;
}

NMOD_INLINE
mp_limb_t nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
   const mp_limb_t neg = mod.n - a;
   if (neg > b)
      return a + b;
   else
      return b - neg;
}

NMOD_INLINE
mp_limb_t nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
   const mp_limb_t diff = a - b;

   if (a < b)
      return mod.n + diff;
   else
      return diff;
}

NMOD_INLINE
mp_limb_t nmod_neg(mp_limb_t a, nmod_t mod)
{
   if (a)
      return mod.n - a;
   else
      return 0;
}

NMOD_INLINE
mp_limb_t nmod_mul(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
    mp_limb_t res;
    NMOD_MUL_PRENORM(res, a, b << mod.norm, mod);
    return res;
}

NMOD_INLINE
mp_limb_t _nmod_mul_fullword(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
    mp_limb_t res;
    NMOD_MUL_FULLWORD(res, a, b, mod);
    return res;
}

NMOD_INLINE
mp_limb_t nmod_addmul(mp_limb_t a, mp_limb_t b, mp_limb_t c, nmod_t mod)
{
    return nmod_add(a, nmod_mul(b, c, mod), mod);
}

#define NMOD_ADDMUL(r, a, b, mod) \
    do { \
       (r) = nmod_addmul((r), (a), (b), (mod)); \
    } while (0)

NMOD_INLINE
mp_limb_t nmod_inv(mp_limb_t a, nmod_t mod)
{
    return n_invmod(a, mod.n);
}

NMOD_INLINE
mp_limb_t nmod_div(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
    return nmod_mul(a, n_invmod(b, mod.n), mod);
}

int nmod_divides(mp_limb_t * a, mp_limb_t b, mp_limb_t c, nmod_t mod);

NMOD_INLINE
mp_limb_t nmod_pow_ui(mp_limb_t a, ulong exp, nmod_t mod)
{
    return n_powmod2_ui_preinv(a, exp, mod.n, mod.ninv);
}

NMOD_INLINE
mp_limb_t nmod_pow_fmpz(mp_limb_t a, const fmpz_t exp, nmod_t mod)
{
    return n_powmod2_fmpz_preinv(a, exp, mod.n, mod.ninv);
}

NMOD_INLINE
void nmod_init(nmod_t * mod, mp_limb_t n)
{
   mod->n = n;
   mod->ninv = n_preinvert_limb(n);
   mod->norm = flint_clz(n);
}

/* discrete logs a la Pohlig - Hellman ***************************************/

typedef struct {
    mp_limb_t gammapow;
    ulong cm;
} nmod_discrete_log_pohlig_hellman_table_entry_struct;

typedef struct {
    slong exp;
    ulong prime;
    mp_limb_t gamma;
    mp_limb_t gammainv;
    mp_limb_t startingbeta;
    ulong co;
    ulong startinge;
    ulong idem;
    ulong cbound;
    ulong dbound;
    nmod_discrete_log_pohlig_hellman_table_entry_struct * table; /* length cbound */
} nmod_discrete_log_pohlig_hellman_entry_struct;

typedef struct {
    nmod_t mod;         /* p is mod.n */
    mp_limb_t alpha;    /* p.r. of p */
    mp_limb_t alphainv;
    slong num_factors;  /* factors of p - 1*/
    nmod_discrete_log_pohlig_hellman_entry_struct * entries;
} nmod_discrete_log_pohlig_hellman_struct;

typedef nmod_discrete_log_pohlig_hellman_struct nmod_discrete_log_pohlig_hellman_t[1];

void nmod_discrete_log_pohlig_hellman_init(
                nmod_discrete_log_pohlig_hellman_t L);

void nmod_discrete_log_pohlig_hellman_clear(
                nmod_discrete_log_pohlig_hellman_t L);

double nmod_discrete_log_pohlig_hellman_precompute_prime(
                nmod_discrete_log_pohlig_hellman_t L,
                mp_limb_t p);

ulong nmod_discrete_log_pohlig_hellman_run(
                const nmod_discrete_log_pohlig_hellman_t L,
                mp_limb_t y);

NMOD_INLINE mp_limb_t nmod_discrete_log_pohlig_hellman_primitive_root(
                const nmod_discrete_log_pohlig_hellman_t L)
{
    return L->alpha;
}

#ifdef __cplusplus
}
#endif

#endif

