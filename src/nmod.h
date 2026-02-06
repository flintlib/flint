/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
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
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NMOD_RED2(r, a_hi, a_lo, mod) \
  do { \
    ulong q0xx, q1xx, r1xx; \
    const ulong u1xx = ((a_hi)<<(mod).norm) \
     + (((mod).norm == 0) ? UWORD(0) : (a_lo)>>(FLINT_BITS - (mod).norm)); \
    const ulong u0xx = (a_lo)<<(mod).norm; \
    const ulong nxx = (mod).n<<(mod).norm; \
    umul_ppmm(q1xx, q0xx, (mod).ninv, u1xx); \
    add_ssaaaa(q1xx, q0xx, q1xx, q0xx, u1xx, u0xx); \
    r1xx = (u0xx - (q1xx + 1)*nxx); \
    if (r1xx > q0xx) r1xx += nxx; \
    if (r1xx < nxx) r = (r1xx>>(mod).norm); \
    else r = ((r1xx - nxx)>>(mod).norm); \
  } while (0)

#define NMOD_RED(r, a, mod) \
  do { \
    NMOD_RED2(r, UWORD(0), a, mod); \
  } while (0)

#define NMOD2_RED2(r, a_hi, a_lo, mod) \
  do { \
    ulong v_hi;	\
    NMOD_RED(v_hi, a_hi, mod); \
    NMOD_RED2(r, v_hi, a_lo, mod); \
  } while (0)

#define NMOD_RED3(r, a_hi, a_me, a_lo, mod) \
  do { \
    ulong v_hi;	\
    NMOD_RED2(v_hi, a_hi, a_me, mod); \
    NMOD_RED2(r, v_hi, a_lo, mod); \
  } while (0)

#define NMOD_BITS(mod) (FLINT_BITS - ((mod).norm))
#define NMOD_CAN_USE_SHOUP(mod) ((mod).norm > 0)

#define NMOD_RED2_FULLWORD(r, a_hi, a_lo, mod) \
  do { \
    ulong q0xx, q1xx, r1xx; \
    const ulong u1xx = (a_hi); \
    const ulong u0xx = (a_lo); \
    const ulong nxx = (mod).n; \
    umul_ppmm(q1xx, q0xx, (mod).ninv, u1xx); \
    add_ssaaaa(q1xx, q0xx, q1xx, q0xx, u1xx, u0xx); \
    r1xx = (u0xx - (q1xx + 1)*nxx); \
    if (r1xx > q0xx) r1xx += nxx; \
    if (r1xx < nxx) r = (r1xx); \
    else r = ((r1xx - nxx)); \
  } while (0)

#define NMOD_RED2_NONFULLWORD(r, a_hi, a_lo, mod) \
  do { \
    ulong q0xx, q1xx, r1xx; \
    const ulong u1xx = ((a_hi)<<(mod).norm) + ((a_lo)>>(FLINT_BITS - (mod).norm)); \
    const ulong u0xx = (a_lo)<<(mod).norm; \
    const ulong nxx = (mod).n<<(mod).norm; \
    umul_ppmm(q1xx, q0xx, (mod).ninv, u1xx); \
    add_ssaaaa(q1xx, q0xx, q1xx, q0xx, u1xx, u0xx); \
    r1xx = (u0xx - (q1xx + 1)*nxx); \
    if (r1xx > q0xx) r1xx += nxx; \
    if (r1xx < nxx) r = (r1xx>>(mod).norm); \
    else r = ((r1xx - nxx)>>(mod).norm); \
  } while (0)

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

NMOD_INLINE ulong nmod_set_ui(ulong x, nmod_t mod)
{
    if (x < mod.n)
        return x;

    NMOD_RED(x, x, mod);
    return x;
}

NMOD_INLINE
ulong nmod_set_si(slong x, nmod_t mod)
{
    ulong res = (x >= 0) ? (ulong) x : -(ulong) x;
    NMOD_RED(res, res, mod);
    return (res == 0 || x > 0) ? res : mod.n - res;
}

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
    FLINT_ASSERT(a < mod.n);
    FLINT_ASSERT(b < mod.n);
    if (neg > b)
        return a + b;
    else
        return b - neg;
}

NMOD_INLINE
ulong nmod_ui_add_ui(ulong a, ulong b, nmod_t mod)
{
    return nmod_add(nmod_set_ui(a, mod), nmod_set_ui(b, mod), mod);
}

NMOD_INLINE
ulong nmod_sub(ulong a, ulong b, nmod_t mod)
{
    const ulong diff = a - b;
    FLINT_ASSERT(a < mod.n);
    FLINT_ASSERT(b < mod.n);
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
    FLINT_ASSERT(a < mod.n);
    FLINT_ASSERT(b < mod.n);
    NMOD_MUL_PRENORM(res, a, b << mod.norm, mod);
    return res;
}

NMOD_INLINE ulong
nmod_ui_mul_ui(ulong a, ulong b, nmod_t mod)
{
    return n_mulmod2_preinv(a, b, mod.n, mod.ninv);
}

NMOD_INLINE
ulong _nmod_mul_fullword(ulong a, ulong b, nmod_t mod)
{
    ulong res;
    NMOD_MUL_FULLWORD(res, a, b, mod);
    return res;
}

NMOD_INLINE
ulong nmod_addmul(ulong s, ulong a, ulong b, nmod_t mod)
{
    ulong hi, lo;
    FLINT_ASSERT(s < mod.n);
    FLINT_ASSERT(a < mod.n);
    FLINT_ASSERT(b < mod.n);
    umul_ppmm(hi, lo, a, b);
    add_ssaaaa(hi, lo, hi, lo, 0, s);
    NMOD_RED2(lo, hi, lo, mod);
    return lo;
}

#define NMOD_ADDMUL(r, a, b, mod) \
    do { \
       (r) = nmod_addmul((r), (a), (b), (mod)); \
    } while (0)

// TODO doc  a*b + c*d
NMOD_INLINE
ulong nmod_fmma(ulong a, ulong b, ulong c, ulong d, nmod_t mod)
{
    a = nmod_mul(a, b, mod);
    NMOD_ADDMUL(a, c, d, mod);
    return a;
}

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

int nmod_divides(ulong * a, ulong b, ulong c, nmod_t mod);

ulong _nmod_pow_ui_redc(ulong a, ulong exp, nmod_t mod);
ulong _nmod_pow_ui_binexp(ulong a, ulong exp, nmod_t mod);
ulong nmod_pow_ui(ulong a, ulong exp, nmod_t mod);

NMOD_INLINE ulong
nmod_ui_pow_ui(ulong a, ulong exp, nmod_t mod)
{
    return nmod_pow_ui(nmod_set_ui(a, mod), exp, mod);
}

ulong _nmod_2_pow_ui_binexp(ulong exp, nmod_t mod);
ulong nmod_2_pow_ui(ulong exp, nmod_t mod);

NMOD_INLINE
ulong nmod_pow_fmpz(ulong a, const fmpz_t exp, nmod_t mod)
{
    return n_powmod2_fmpz_preinv(a, exp, mod.n, mod.ninv);
}

NMOD_INLINE
void nmod_init(nmod_t * mod, ulong n)
{
    mod->n = n;
    mod->norm = flint_clz(n);
    mod->ninv = n_preinvert_limb_prenorm(n << (mod->norm));
}

/* Montgomery arithmetic *****************************************************/

/* Internal helpers for Montgomery arithmetic ********************************/
/* Some of these functions could be moved to ulong_extras in the future. */

/* Experimental: see if the compiler generates better code than inline assembly */
#if FLINT_BITS == 64 && defined(__GNUC__)

typedef __uint128_t ull_t;

NMOD_INLINE ulong ull_hi(ull_t x) { return x >> FLINT_BITS; }
NMOD_INLINE ulong ull_lo(ull_t x) { return (ulong) x; }
NMOD_INLINE ull_t ull_add(ull_t x, ull_t y) { return x + y; }
NMOD_INLINE ull_t ull_add_u(ull_t x, ulong y) { return x + (ull_t) y; }
NMOD_INLINE ull_t ull_u_mul_u(ulong x, ulong y) { return (ull_t) x * (ull_t) y; }
NMOD_INLINE ull_t ull(ulong hi, ulong lo) { return ((ull_t) lo) | (((ull_t) hi) << FLINT_BITS); }

#else

typedef struct { ulong lo; ulong hi; } ull_t;

NMOD_INLINE ulong ull_hi(ull_t x) { return x.hi; }
NMOD_INLINE ulong ull_lo(ull_t x) { return x.lo; }
NMOD_INLINE ull_t ull_add(ull_t x, ull_t y) { ull_t t; add_ssaaaa(t.hi, t.lo, x.hi, x.lo, y.hi, y.lo); return t; }
NMOD_INLINE ull_t ull_add_u(ull_t x, ulong y) { ull_t t; add_ssaaaa(t.hi, t.lo, x.hi, x.lo, 0, y); return t; }
NMOD_INLINE ull_t ull_u_mul_u(ulong x, ulong y) { ull_t t; umul_ppmm(t.hi, t.lo, x, y); return t; }
NMOD_INLINE ull_t ull(ulong hi, ulong lo) { ull_t t; t.lo = lo; t.hi = hi; return t; }

#endif

FLINT_FORCE_INLINE ulong n_mulhi(ulong a, ulong b)
{
    return ull_hi(ull_u_mul_u(a, b));
}

FLINT_FORCE_INLINE void n_mul2(ulong * hi, ulong * lo, ulong a, ulong b)
{
    ull_t t = ull_u_mul_u(a, b);
    *hi = ull_hi(t);
    *lo = ull_lo(t);
}

/* Assumes a already reduced mod n. */
FLINT_FORCE_INLINE ulong
n_to_redc_preinv(ulong a, ulong n, ulong ninv, ulong norm)
{
    ulong q0, q1, r;
    FLINT_ASSERT(a < n);
    n <<= norm;
    a <<= norm;
    ull_t t = ull_u_mul_u(ninv, a);
    q1 = ull_hi(t);
    q0 = ull_lo(t);
    q1 += a;
    r = -(q1 + 1) * n;
    if (r > q0)
        r += n;
    return (r < n) ? (r >> norm) : ((r - n) >> norm);
}

/* Computes x/R, reducing to [0,2n). */
/* nred = -1/n mod R, R = 2^FLINT_BITS */
FLINT_FORCE_INLINE
ulong n_redc_fast(ulong x, ulong n, ulong nred)
{
    return ull_hi(ull_add_u(ull_u_mul_u(n, x * nred), x));
}

FLINT_FORCE_INLINE
ulong n_ll_redc_fast(ull_t x, ulong n, ulong nred)
{
    return ull_hi(ull_add(ull_u_mul_u(n, ull_lo(x) * nred), x));
}

/* Computes x/R mod n, reducing to [0,n). */
FLINT_FORCE_INLINE
ulong n_redc(ulong x, ulong n, ulong nred)
{
    ulong y = n_redc_fast(x, n, nred);
    if (y >= n)
        y -= n;
    return y;
}

/* Computes x/R, reducing to [0,n). Assumes x < n * 2^FLINT_BITS. */
FLINT_FORCE_INLINE
ulong n_ll_redc(ull_t x, ulong n, ulong nred)
{
    ulong lo, hi;

    lo = n_mulhi(ull_lo(x) * (-nred), n);
    hi = ull_hi(x);

    if (hi < lo)
        return hi - lo + n;
    else
        return hi - lo;
}

/* Computes (ab)/R mod n, reducing to [0,n). */
FLINT_FORCE_INLINE
ulong n_mulmod_redc(ulong a, ulong b, ulong n, ulong nred)
{
    return n_ll_redc(ull_u_mul_u(a, b), n, nred);
}

/* Computes (ab)/R mod n, reducing to [0,2n). */
FLINT_FORCE_INLINE
ulong n_mulmod_redc_fast(ulong a, ulong b, ulong n, ulong nred)
{
    return n_ll_redc_fast(ull_u_mul_u(a, b), n, nred);
}

/* User-friendly functions for Montgomery arithmetic *************************/

typedef struct
{
    /* Standard nmod data for odd n. The precomputed inverse is not required
       for Montgomery arithmetic itself, but it is useful for conversion
       into Montgomery form. */
    /* nred = -1/n mod R, R = 2^FLINT_BITS */
    nmod_t mod;
    ulong nred;
    /* Todo: should we store the constant 1? This would slow down
       creating the context when we don't need it.
    ulong one;
    */
}
nmod_redc_ctx_struct;

typedef nmod_redc_ctx_struct nmod_redc_ctx_t[1];

NMOD_INLINE void nmod_redc_ctx_init_nmod(nmod_redc_ctx_t ctx, nmod_t mod)
{
    FLINT_ASSERT(mod.n & 1);

    ctx->mod = mod;
    ctx->nred = n_binvert(-mod.n);
}

NMOD_INLINE void nmod_redc_ctx_init_ui(nmod_redc_ctx_t ctx, ulong n)
{
    FLINT_ASSERT(n & 1);

    nmod_init(&ctx->mod, n);
    ctx->nred = n_binvert(-n);
}

NMOD_INLINE ulong nmod_redc_set_nmod(ulong x, const nmod_redc_ctx_t ctx)
{
    return n_to_redc_preinv(x, ctx->mod.n, ctx->mod.ninv, ctx->mod.norm);
}

NMOD_INLINE ulong nmod_redc_set_ui(ulong x, const nmod_redc_ctx_t ctx)
{
    return nmod_redc_set_nmod(nmod_set_ui(x, ctx->mod), ctx);
}

NMOD_INLINE ulong nmod_redc_get_nmod(ulong x, const nmod_redc_ctx_t ctx)
{
    return n_redc(x, ctx->mod.n, ctx->nred);
}

NMOD_INLINE ulong nmod_redc_neg(ulong x, const nmod_redc_ctx_t ctx)
{
    return nmod_neg(x, ctx->mod);
}

NMOD_INLINE ulong nmod_redc_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return nmod_add(x, y, ctx->mod);
}

NMOD_INLINE ulong nmod_redc_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return nmod_sub(x, y, ctx->mod);
}

NMOD_INLINE ulong nmod_redc_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return n_mulmod_redc(x, y, ctx->mod.n, ctx->nred);
}

NMOD_INLINE int nmod_redc_can_use_fast(const nmod_redc_ctx_t ctx)
{
    return ctx->mod.norm >= 2;
}

NMOD_INLINE ulong nmod_redc_fast_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return n_mulmod_redc_fast(x, y, ctx->mod.n, ctx->nred);
}

NMOD_INLINE ulong nmod_redc_fast_neg(ulong x, const nmod_redc_ctx_t ctx)
{
    return n_negmod(x, 2 * ctx->mod.n);
}

NMOD_INLINE ulong nmod_redc_fast_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return n_addmod(x, y, 2 * ctx->mod.n);
}

NMOD_INLINE ulong nmod_redc_fast_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return n_submod(x, y, 2 * ctx->mod.n);
}

NMOD_INLINE ulong nmod_redc_fast_normalise(ulong x, const nmod_redc_ctx_t ctx)
{
    return x < ctx->mod.n ? x : x - ctx->mod.n;
}

NMOD_INLINE
ulong nmod_redc_fast_mul_two(ulong x, const nmod_redc_ctx_t ctx)
{
    ulong n = ctx->mod.n;

    if (x >= n)
        x -= n;

    return 2 * x;
}

ulong _nmod_redc_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx);
ulong _nmod_redc_fast_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx);

ulong _nmod_redc_2_pow_ui(ulong exp, const nmod_redc_ctx_t ctx);
ulong _nmod_redc_fast_2_pow_ui(ulong exp, const nmod_redc_ctx_t ctx);

/* Half-limb versions of Montgomery arithmetic */

FLINT_FORCE_INLINE ulong
n_to_redc_half_preinv(ulong a, ulong n, ulong ninv, ulong FLINT_UNUSED(norm))
{
    return n_mod2_preinv(((ulong) a) << (FLINT_BITS / 2), n, ninv);
}

FLINT_FORCE_INLINE
ulong n_redc_half_fast(ulong x, ulong n, ulong nred)
{
    ulong y = (x * nred) % (UWORD(1) << (FLINT_BITS / 2));
    ulong z = x + n * y;
    return z >> (FLINT_BITS / 2);
}

FLINT_FORCE_INLINE
ulong n_redc_half(ulong x, ulong n, ulong nred)
{
    ulong y = n_redc_half_fast(x, n, nred);
    if (y >= n)
        y -= n;
    return y;
}

FLINT_FORCE_INLINE
ulong n_mulmod_redc_half(ulong a, ulong b, ulong n, ulong nred)
{
    return n_redc_half(a * b, n, nred);
}

FLINT_FORCE_INLINE
ulong n_mulmod_redc_half_fast(ulong a, ulong b, ulong n, ulong nred)
{
    return n_redc_half_fast(a * b, n, nred);
}

NMOD_INLINE void nmod_redc_half_ctx_init_nmod(nmod_redc_ctx_t ctx, nmod_t mod)
{
    FLINT_ASSERT(mod.n & 1);
    FLINT_ASSERT(mod.norm >= FLINT_BITS / 2 + 1);
    ctx->mod = mod;
    /* todo: skip one iteration */
    ctx->nred = n_binvert(-mod.n) % (UWORD(1) << (FLINT_BITS / 2));
}

NMOD_INLINE void nmod_redc_half_ctx_init_ui(nmod_redc_ctx_t ctx, ulong n)
{
    nmod_init(&ctx->mod, n);
    nmod_redc_half_ctx_init_nmod(ctx, ctx->mod);
}

NMOD_INLINE ulong nmod_redc_half_set_nmod(ulong x, const nmod_redc_ctx_t ctx)
{
    return n_to_redc_half_preinv(x, ctx->mod.n, ctx->mod.ninv, ctx->mod.norm);
}

NMOD_INLINE ulong nmod_redc_half_set_ui(ulong x, const nmod_redc_ctx_t ctx)
{
    return nmod_redc_half_set_nmod(nmod_set_ui(x, ctx->mod), ctx);
}

NMOD_INLINE ulong nmod_redc_half_get_nmod(ulong x, const nmod_redc_ctx_t ctx)
{
    return n_redc_half(x, ctx->mod.n, ctx->nred);
}

NMOD_INLINE ulong nmod_redc_half_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return nmod_add(x, y, ctx->mod);
}

NMOD_INLINE ulong nmod_redc_half_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return nmod_sub(x, y, ctx->mod);
}

NMOD_INLINE ulong nmod_redc_half_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return n_mulmod_redc_half(x, y, ctx->mod.n, ctx->nred);
}

NMOD_INLINE int nmod_redc_half_can_use_fast(const nmod_redc_ctx_t ctx)
{
    return ctx->mod.norm >= (FLINT_BITS / 2 + 2);
}

NMOD_INLINE ulong nmod_redc_half_fast_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return n_mulmod_redc_half_fast(x, y, ctx->mod.n, ctx->nred);
}

NMOD_INLINE ulong nmod_redc_half_fast_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return n_addmod(x, y, 2 * ctx->mod.n);
}

NMOD_INLINE ulong nmod_redc_half_fast_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
{
    return n_submod(x, y, 2 * ctx->mod.n);
}

/* Rings */

typedef struct
{
    nmod_t nmod;
    ulong a;   /* when used as finite field with defining polynomial x - a */
    truth_t is_prime;
}
_gr_nmod_ctx_struct;

#define NMOD_CTX_REF(ring_ctx) (&((((_gr_nmod_ctx_struct *)(ring_ctx))->nmod)))
#define NMOD_CTX(ring_ctx) (*NMOD_CTX_REF(ring_ctx))
#define NMOD_IS_PRIME(ring_ctx) (((_gr_nmod_ctx_struct *)(ring_ctx))->is_prime)
/* when used as finite field when defining polynomial x - a, allow storing the coefficient a */
#define NMOD_CTX_A(ring_ctx) (&((((_gr_nmod_ctx_struct *)(ring_ctx))->a)))


typedef struct
{
    nmod_redc_ctx_struct ctx;
    ulong one;
    truth_t is_prime;
}
_gr_nmod_redc_ctx_struct;

#define GR_NMOD_REDC_CTX(ring_ctx) (&((((_gr_nmod_redc_ctx_struct *)(ring_ctx))->ctx)))
#define GR_NMOD_REDC_MOD(ring_ctx) (((((_gr_nmod_redc_ctx_struct *)(ring_ctx))->ctx)).mod)
#define GR_NMOD_REDC_N(ring_ctx) ((((((_gr_nmod_redc_ctx_struct *)(ring_ctx))->ctx)).mod).n)
#define GR_NMOD_REDC_NRED(ring_ctx) (GR_NMOD_REDC_CTX(ring_ctx)->nred)
#define GR_NMOD_REDC_ONE(ring_ctx) ((((_gr_nmod_redc_ctx_struct *)(ring_ctx))->one))
#define GR_NMOD_REDC_IS_PRIME(ring_ctx) ((((_gr_nmod_redc_ctx_struct *)(ring_ctx))->is_prime))
#define GR_NMOD_REDC_IS_FAST(ring_ctx) ((ring_ctx)->which_ring == GR_CTX_NMOD_REDC_FAST)

int gr_ctx_init_nmod_redc(gr_ctx_t ctx, ulong n);
int gr_ctx_init_nmod_redc_fast(gr_ctx_t ctx, ulong n);

/* Discrete logs a la Pohlig - Hellman ***************************************/

typedef struct {
    ulong gammapow;
    ulong cm;
} nmod_discrete_log_pohlig_hellman_table_entry_struct;

typedef struct {
    slong exp;
    ulong prime;
    ulong gamma;
    ulong gammainv;
    ulong startingbeta;
    ulong co;
    ulong startinge;
    ulong idem;
    ulong cbound;
    ulong dbound;
    nmod_discrete_log_pohlig_hellman_table_entry_struct * table; /* length cbound */
} nmod_discrete_log_pohlig_hellman_entry_struct;

typedef struct {
    nmod_t mod;         /* p is mod.n */
    ulong alpha;    /* p.r. of p */
    ulong alphainv;
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
                ulong p);

ulong nmod_discrete_log_pohlig_hellman_run(
                const nmod_discrete_log_pohlig_hellman_t L,
                ulong y);

NMOD_INLINE ulong nmod_discrete_log_pohlig_hellman_primitive_root(
                const nmod_discrete_log_pohlig_hellman_t L)
{
    return L->alpha;
}

#ifdef __cplusplus
}
#endif

#endif
