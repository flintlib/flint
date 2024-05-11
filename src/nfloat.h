/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NFLOAT_H
#define NFLOAT_H

#ifdef NFLOAT_INLINES_C
#define NFLOAT_INLINE
#else
#define NFLOAT_INLINE static inline
#endif

#include "flint.h"
#include "mpn_extras.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_mat.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The upper limit is so that temporary buffers are safe to allocate on
   the stack and so that simple operations like swapping are not too
   expensive compared to a pointer-and-size representation. */
#if FLINT_BITS == 64
#define NFLOAT_MIN_LIMBS 1
#define NFLOAT_MAX_LIMBS 33
#else
#define NFLOAT_MIN_LIMBS 1
#define NFLOAT_MAX_LIMBS 66
#endif

/* Number of header limbs used to encode sign + exponent. We use a
   whole limb for the sign bit to avoid the overhead of bit fiddling. */
#define NFLOAT_HEADER_LIMBS 2
#define NFLOAT_EXP(x) (((mp_limb_signed_t *) x)[0])
#define NFLOAT_SGNBIT(x) (((mp_ptr) x)[1])
#define NFLOAT_D(x) (((mp_ptr) x) + NFLOAT_HEADER_LIMBS)
#define NFLOAT_DATA(x) (((mp_ptr) x))

/* Limbs needed to hold a temporary element of any precision. */
#define NFLOAT_MAX_ALLOC (NFLOAT_HEADER_LIMBS + NFLOAT_MAX_LIMBS)

/* The exponent limits should not be too close to the word boundary,
   so that we can add two exponents plus some adjustmenets
   and safely intercept overflow/underflow. */
#define NFLOAT_MIN_EXP (WORD_MIN / 4)
#define NFLOAT_MAX_EXP (WORD_MAX / 4)

/* Special values (including zero). */
/* TODO: is it better to encode special values in the exponent or in
         the sign field? */
/* TODO: do we want signed zero (optionally)? */
#define NFLOAT_EXP_ZERO      WORD_MIN
#define NFLOAT_EXP_POS_INF   (WORD_MIN+1)
#define NFLOAT_EXP_NEG_INF   (WORD_MIN+2)
#define NFLOAT_EXP_NAN       (WORD_MIN+3)
#define NFLOAT_IS_SPECIAL(x) (NFLOAT_EXP(x) < NFLOAT_MIN_EXP)
#define NFLOAT_IS_ZERO(x)    (NFLOAT_EXP(x) == NFLOAT_EXP_ZERO)
#define NFLOAT_IS_POS_INF(x) (NFLOAT_EXP(x) == NFLOAT_EXP_POS_INF)
#define NFLOAT_IS_NEG_INF(x) (NFLOAT_EXP(x) == NFLOAT_EXP_NEG_INF)
#define NFLOAT_IS_INF(x)     (NFLOAT_IS_POS_INF(x) || NFLOAT_IS_NEG_INF(x))
#define NFLOAT_IS_NAN(x)     (NFLOAT_EXP(x) == NFLOAT_EXP_NAN)

/* Context flags to indicate whether we want to allow special values
   or return GR_DOMAIN / GR_UNABLE. */
#define NFLOAT_ALLOW_UNDERFLOW 1
#define NFLOAT_ALLOW_INF       2
#define NFLOAT_ALLOW_NAN       4

typedef struct
{
    slong nlimbs;
    int flags;
    int rnd;    /* Allow rounding modes? Currently unused. */
}
_nfloat_ctx_struct;

typedef void * nfloat_ptr;
typedef const void * nfloat_srcptr;

#define NFLOAT_CTX(ctx) ((_nfloat_ctx_struct *)(ctx))
#define NFLOAT_CTX_NLIMBS(ctx) (NFLOAT_CTX(ctx)->nlimbs)
#define NFLOAT_CTX_PREC(ctx) ((NFLOAT_CTX(ctx)->nlimbs) * FLINT_BITS)
#define NFLOAT_CTX_FLAGS(ctx) (NFLOAT_CTX(ctx)->flags)
#define NFLOAT_CTX_RND(ctx) (NFLOAT_CTX(ctx)->rnd)
#define NFLOAT_CTX_DATA_NLIMBS(ctx) (NFLOAT_CTX_NLIMBS(ctx) + NFLOAT_HEADER_LIMBS)
#define NFLOAT_CTX_HAS_INF_NAN(ctx) ((NFLOAT_CTX_FLAGS(ctx) & (NFLOAT_ALLOW_INF | NFLOAT_ALLOW_NAN)) != 0)

typedef struct { mp_limb_t head[NFLOAT_HEADER_LIMBS]; mp_limb_t d[64 / FLINT_BITS]; } nfloat64_struct;
typedef struct { mp_limb_t head[NFLOAT_HEADER_LIMBS]; mp_limb_t d[128 / FLINT_BITS]; } nfloat128_struct;
typedef struct { mp_limb_t head[NFLOAT_HEADER_LIMBS]; mp_limb_t d[192 / FLINT_BITS]; } nfloat192_struct;
typedef struct { mp_limb_t head[NFLOAT_HEADER_LIMBS]; mp_limb_t d[256 / FLINT_BITS]; } nfloat256_struct;
typedef struct { mp_limb_t head[NFLOAT_HEADER_LIMBS]; mp_limb_t d[384 / FLINT_BITS]; } nfloat384_struct;
typedef struct { mp_limb_t head[NFLOAT_HEADER_LIMBS]; mp_limb_t d[512 / FLINT_BITS]; } nfloat512_struct;
typedef struct { mp_limb_t head[NFLOAT_HEADER_LIMBS]; mp_limb_t d[1024 / FLINT_BITS]; } nfloat1024_struct;
typedef struct { mp_limb_t head[NFLOAT_HEADER_LIMBS]; mp_limb_t d[2048 / FLINT_BITS]; } nfloat2048_struct;

typedef nfloat64_struct nfloat64_t[1];
typedef nfloat128_struct nfloat128_t[1];
typedef nfloat192_struct nfloat192_t[1];
typedef nfloat256_struct nfloat256_t[1];
typedef nfloat384_struct nfloat384_t[1];
typedef nfloat512_struct nfloat512_t[1];
typedef nfloat1024_struct nfloat1024_t[1];
typedef nfloat2048_struct nfloat2048_t[1];

#define LIMB_MSB_IS_SET(n) ((mp_limb_signed_t) (n) < 0)

int nfloat_ctx_init(gr_ctx_t ctx, slong prec, int flags);
int nfloat_ctx_write(gr_stream_t out, gr_ctx_t ctx);

NFLOAT_INLINE void
nfloat_init(nfloat_ptr res, gr_ctx_t ctx)
{
    NFLOAT_EXP(res) = NFLOAT_EXP_ZERO;
    NFLOAT_SGNBIT(res) = 0;
}

NFLOAT_INLINE void
nfloat_clear(nfloat_ptr res, gr_ctx_t ctx)
{
}

void nfloat_swap(nfloat_ptr x, nfloat_ptr y, gr_ctx_t ctx);
int nfloat_set(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);

truth_t nfloat_equal(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx);

NFLOAT_INLINE
int _nfloat_ctx_set_real_prec(gr_ctx_t ctx, slong prec)
{
    return GR_UNABLE;
}

NFLOAT_INLINE
int _nfloat_ctx_get_real_prec(slong * res, gr_ctx_t ctx)
{
    *res = NFLOAT_CTX_PREC(ctx);
    return GR_SUCCESS;
}

NFLOAT_INLINE int
nfloat_zero(nfloat_ptr res, gr_ctx_t ctx)
{
    NFLOAT_EXP(res) = NFLOAT_EXP_ZERO;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}

int nfloat_pos_inf(nfloat_ptr res, gr_ctx_t ctx);
int nfloat_neg_inf(nfloat_ptr res, gr_ctx_t ctx);
int nfloat_nan(nfloat_ptr res, gr_ctx_t ctx);

int _nfloat_underflow(nfloat_ptr res, int sgnbit, gr_ctx_t ctx);
int _nfloat_overflow(nfloat_ptr res, int sgnbit, gr_ctx_t ctx);

#define NFLOAT_HANDLE_UNDERFLOW(res, ctx) \
    do { \
        if (FLINT_UNLIKELY(NFLOAT_EXP(res) < NFLOAT_MIN_EXP)) \
            return _nfloat_underflow(res, NFLOAT_SGNBIT(res), ctx); \
    } while (0)

#define NFLOAT_HANDLE_OVERFLOW(res, ctx) \
    do { \
        if (FLINT_UNLIKELY(NFLOAT_EXP(res) < NFLOAT_MIN_EXP)) \
            return _nfloat_underflow(res, NFLOAT_SGNBIT(res), ctx); \
    } while (0)

#define NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx) \
    do { \
        NFLOAT_HANDLE_UNDERFLOW(res, ctx); \
        NFLOAT_HANDLE_OVERFLOW(res, ctx); \
    } while (0)

int nfloat_one(nfloat_ptr res, gr_ctx_t ctx);
int nfloat_neg_one(nfloat_ptr res, gr_ctx_t ctx);

NFLOAT_INLINE truth_t
nfloat_is_zero(nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_NAN(x))
        return T_UNKNOWN;

    return NFLOAT_IS_ZERO(x) ? T_TRUE : T_FALSE;
}

truth_t nfloat_is_one(nfloat_srcptr x, gr_ctx_t ctx);
truth_t nfloat_is_neg_one(nfloat_srcptr x, gr_ctx_t ctx);

int nfloat_set_ui(nfloat_ptr res, ulong x, gr_ctx_t ctx);
int nfloat_set_si(nfloat_ptr res, slong x, gr_ctx_t ctx);

/* Here exp is understood such that {x, xn} is a fraction in [0, 1). */
NFLOAT_INLINE int
_nfloat_set_mpn_2exp(nfloat_ptr res, mp_srcptr x, mp_size_t xn, slong exp, int xsgnbit, gr_ctx_t ctx)
{
    mp_limb_t top;
    slong norm;
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(x[xn - 1] != 0);

    top = x[xn - 1];

    if (LIMB_MSB_IS_SET(top))
    {
        if (xn >= nlimbs)
        {
            flint_mpn_copyi(NFLOAT_D(res), x + xn - nlimbs, nlimbs);
        }
        else
        {
            flint_mpn_zero(NFLOAT_D(res), nlimbs - xn);
            flint_mpn_copyi(NFLOAT_D(res) + nlimbs - xn, x, xn);
        }
    }
    else
    {
        norm = flint_clz(top);

        if (xn > nlimbs)
        {
            mpn_lshift(NFLOAT_D(res), x + xn - nlimbs, nlimbs, norm);
            NFLOAT_D(res)[0] |= (x[xn - nlimbs - 1] >> (FLINT_BITS - norm));
        }
        else
        {
            flint_mpn_zero(NFLOAT_D(res), nlimbs - xn);
            mpn_lshift(NFLOAT_D(res) + nlimbs - xn, x, xn, norm);
        }

        exp -= norm;
    }

    NFLOAT_SGNBIT(res) = xsgnbit;
    NFLOAT_EXP(res) = exp;
    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);

    return GR_SUCCESS;
}

NFLOAT_INLINE int
nfloat_set_mpn_2exp(nfloat_ptr res, mp_srcptr x, mp_size_t xn, slong exp, int xsgnbit, gr_ctx_t ctx)
{
    while (xn != 0 && x[xn - 1] == 0)
    {
        xn--;
        exp -= FLINT_BITS;
    }

    if (xn == 0)
        return nfloat_zero(res, ctx);

    return _nfloat_set_mpn_2exp(res, x, xn, exp, xsgnbit, ctx);
}

int nfloat_set_fmpz(nfloat_ptr res, const fmpz_t x, gr_ctx_t ctx);

#ifdef ARF_H
int nfloat_set_arf(nfloat_ptr res, const arf_t x, gr_ctx_t ctx);
int nfloat_get_arf(arf_t res, nfloat_srcptr x, gr_ctx_t ctx);
#endif

int nfloat_set_fmpq(nfloat_ptr res, const fmpq_t v, gr_ctx_t ctx);
int nfloat_set_d(nfloat_ptr res, double x, gr_ctx_t ctx);
int nfloat_set_str(nfloat_ptr res, const char * x, gr_ctx_t ctx);
int nfloat_set_other(nfloat_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);

int nfloat_write(gr_stream_t out, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_randtest(nfloat_ptr res, flint_rand_t state, gr_ctx_t ctx);

int nfloat_neg(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_abs(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int _nfloat_cmp(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx);
int _nfloat_cmpabs(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx);
int nfloat_cmp(int * res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx);
int nfloat_cmpabs(int * res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx);
int nfloat_sgn(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_im(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);

int _nfloat_add_1(nfloat_ptr res, mp_limb_t x0, slong xexp, int xsgnbit, mp_limb_t y0, slong delta, gr_ctx_t ctx);
int _nfloat_sub_1(nfloat_ptr res, mp_limb_t x0, slong xexp, int xsgnbit, mp_limb_t y0, slong delta, gr_ctx_t ctx);
int _nfloat_add_n(nfloat_ptr res, mp_srcptr xd, slong xexp, int xsgnbit, mp_srcptr yd, slong delta, slong nlimbs, gr_ctx_t ctx);
int _nfloat_sub_n(nfloat_ptr res, mp_srcptr xd, slong xexp, int xsgnbit, mp_srcptr yd, slong delta, slong nlimbs, gr_ctx_t ctx);

int nfloat_add(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx);
int nfloat_sub(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx);
int nfloat_mul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx);

int nfloat_mul_2exp_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx);

int nfloat_inv(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx);
int nfloat_div_ui(nfloat_ptr res, nfloat_srcptr x, ulong y, gr_ctx_t ctx);
int nfloat_div_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx);

int nfloat_sqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_rsqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);

int nfloat_floor(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_ceil(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_trunc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_nint(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);

int nfloat_pi(nfloat_ptr res, gr_ctx_t ctx);
int nfloat_pow(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx);
int nfloat_exp(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_expm1(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_log(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_log1p(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_sin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_cos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_tan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_sinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_cosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_tanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_atan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_gamma(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);
int nfloat_zeta(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx);

void _nfloat_vec_init(nfloat_ptr res, slong len, gr_ctx_t ctx);
void _nfloat_vec_clear(nfloat_ptr res, slong len, gr_ctx_t ctx);
int _nfloat_vec_set(nfloat_ptr res, nfloat_srcptr x, slong len, gr_ctx_t ctx);
int _nfloat_vec_zero(nfloat_ptr res, slong len, gr_ctx_t ctx);

int _nfloat_vec_add(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx);
int _nfloat_vec_sub(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx);
int _nfloat_vec_mul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx);
int _nfloat_vec_mul_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx);

int _nfloat_vec_dot(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx);
int _nfloat_vec_dot_rev(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
