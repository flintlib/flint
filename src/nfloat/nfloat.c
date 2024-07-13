/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "fmpz.h"
#include "arf.h"
#include "arb.h"
#include "nfloat.h"
#include "gr_generic.h"
#include "gr_special.h"

/* todo: define in longlong.h */
#if FLINT_BITS == 64 && defined(__GNUC__) && defined(__AVX2__)

#define add_sssssaaaaaaaaaa(s4,s3,s2,s1,s0, a4,a3,a2,a1,a0, b4,b3,b2,b1,b0)  \
  __asm__ ("addq %14,%q4\n\tadcq %12,%q3\n\tadcq %10,%q2\n\tadcq %8,%q1\n\tadcq %6,%q0"    \
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((ulong)(a4)), "rme" ((ulong)(b4)),                 \
         "1"  ((ulong)(a3)), "rme" ((ulong)(b3)),                 \
         "2"  ((ulong)(a2)), "rme" ((ulong)(b2)),                 \
         "3"  ((ulong)(a1)), "rme" ((ulong)(b1)),                 \
         "4"  ((ulong)(a0)), "rme" ((ulong)(b0)))


#define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)  \
  __asm__ ("subq %11,%q3\n\tsbbq %9,%q2\n\tsbbq %7,%q1\n\tsbbq %5,%q0"    \
       : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((ulong)(a3)), "rme" ((ulong)(b3)),                 \
         "1"  ((ulong)(a2)), "rme" ((ulong)(b2)),                 \
         "2"  ((ulong)(a1)), "rme" ((ulong)(b1)),                 \
         "3"  ((ulong)(a0)), "rme" ((ulong)(b0)))

#define sub_dddddmmmmmsssss(s4,s3,s2,s1,s0, a4,a3,a2,a1,a0, b4,b3,b2,b1,b0)  \
  __asm__ ("subq %14,%q4\n\tsbbq %12,%q3\n\tsbbq %10,%q2\n\tsbbq %8,%q1\n\tsbbq %6,%q0"    \
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((ulong)(a4)), "rme" ((ulong)(b4)),                 \
         "1"  ((ulong)(a3)), "rme" ((ulong)(b3)),                 \
         "2"  ((ulong)(a2)), "rme" ((ulong)(b2)),                 \
         "3"  ((ulong)(a1)), "rme" ((ulong)(b1)),                 \
         "4"  ((ulong)(a0)), "rme" ((ulong)(b0)))
#else

#define add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0)         \
  do {                                                                                          \
    ulong __t0 = 0;                                                                         \
    add_ssssaaaaaaaa(__t0, s2, s1, s0, (ulong) 0, a2, a1, a0, (ulong) 0, b2, b1, b0);   \
    add_ssaaaa(s4, s3, a4, a3, b4, b3);                                                         \
    add_ssaaaa(s4, s3, s4, s3, (ulong) 0, __t0);                                            \
  } while (0)


#define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)        \
  do {                                                                          \
    ulong __t1, __u1;                                                       \
    sub_dddmmmsss(__t1, s1, s0, (ulong) 0, a1, a0, (ulong) 0, b1, b0);  \
    sub_ddmmss(__u1, s2, (ulong) 0, a2, (ulong) 0, b2);                 \
    sub_ddmmss(s3, s2, (a3) - (b3), s2, -__u1, -__t1);                          \
  } while (0)

#define sub_dddddmmmmmsssss(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0)         \
  do {                                                                                          \
    ulong __t2, __u2;                                                                       \
    sub_ddddmmmmssss(__t2, s2, s1, s0, (ulong) 0, a2, a1, a0, (ulong) 0, b2, b1, b0);   \
    sub_ddmmss(__u2, s3, (ulong) 0, a3, (ulong) 0, b3);                                 \
    sub_ddmmss(s4, s3, (a4) - (b4), s3, -__u2, -__t2);                                          \
  } while (0)

#endif

int
nfloat_write(gr_stream_t out, nfloat_srcptr x, gr_ctx_t ctx)
{
    gr_ctx_t arf_ctx;
    arf_t t;
    int status;

    gr_ctx_init_real_float_arf(arf_ctx, NFLOAT_CTX_PREC(ctx));
    arf_init(t);
    nfloat_get_arf(t, x, ctx);
    status = gr_write(out, t, arf_ctx);
    arf_clear(t);
    return status;
    gr_ctx_clear(arf_ctx);
}

int
nfloat_randtest(nfloat_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    arf_t t;
    int status;

    arf_init(t);
    arf_randtest(t, state, NFLOAT_CTX_PREC(ctx), n_randint(state, 2) ? 2 : 10);
    status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);
    return status;
}

void
nfloat_swap(nfloat_ptr x, nfloat_ptr y, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        FLINT_SWAP(ulong, NFLOAT_DATA(x)[i], NFLOAT_DATA(y)[i]);
}

int
nfloat_set(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        NFLOAT_DATA(res)[i] = NFLOAT_DATA(x)[i];

    return GR_SUCCESS;
}

void
_nfloat_vec_init(nfloat_ptr res, slong len, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < len; i++)
        nfloat_init((nn_ptr) res + i * n, ctx);
}

void
_nfloat_vec_clear(nfloat_ptr res, slong len, gr_ctx_t ctx)
{
}

int
_nfloat_vec_set(nfloat_ptr res, nfloat_srcptr x, slong len, gr_ctx_t ctx)
{
    flint_mpn_copyi(res, x, len * NFLOAT_CTX_DATA_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
_nfloat_vec_zero(nfloat_ptr res, slong len, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < len; i++)
        nfloat_zero((nn_ptr) res + i * n, ctx);

    return GR_SUCCESS;
}

truth_t
nfloat_equal(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        return (NFLOAT_EXP(x) == NFLOAT_EXP(y) &&
                NFLOAT_SGNBIT(x) == NFLOAT_SGNBIT(y)) ? T_TRUE : T_FALSE;
    }

    if (NFLOAT_EXP(x) != NFLOAT_EXP(y))
        return T_FALSE;

    if (NFLOAT_SGNBIT(x) != NFLOAT_SGNBIT(y))
        return T_FALSE;

    for (i = 0; i < n; i++)
        if (NFLOAT_D(x)[i] != NFLOAT_D(y)[i])
            return T_FALSE;

    return T_TRUE;
}

int
nfloat_one(nfloat_ptr res, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    NFLOAT_EXP(res) = 1;
    NFLOAT_SGNBIT(res) = 0;
    flint_mpn_zero(NFLOAT_D(res), n - 1);
    NFLOAT_D(res)[n - 1] = UWORD(1) << (FLINT_BITS - 1);
    return GR_SUCCESS;
}

int
nfloat_neg_one(nfloat_ptr res, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    NFLOAT_EXP(res) = 1;
    NFLOAT_SGNBIT(res) = 1;
    flint_mpn_zero(NFLOAT_D(res), n - 1);
    NFLOAT_D(res)[n - 1] = UWORD(1) << (FLINT_BITS - 1);
    return GR_SUCCESS;
}

truth_t
nfloat_is_one(nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);

    if (NFLOAT_IS_NAN(x))
        return T_UNKNOWN;

    return ((NFLOAT_EXP(x) == 1) && (NFLOAT_SGNBIT(x) == 0) &&
           flint_mpn_zero_p(NFLOAT_D(x), n - 1) &&
            (NFLOAT_D(x)[n - 1] == (UWORD(1) << (FLINT_BITS - 1)))) ? T_TRUE : T_FALSE;
}

truth_t
nfloat_is_neg_one(nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);

    if (NFLOAT_IS_NAN(x))
        return T_UNKNOWN;

    return (NFLOAT_EXP(x) == 1) && (NFLOAT_SGNBIT(x) == 1) &&
           flint_mpn_zero_p(NFLOAT_D(x), n - 1) &&
            (NFLOAT_D(x)[n - 1] == (UWORD(1) << (FLINT_BITS - 1))) ? T_TRUE : T_FALSE;
}

int
nfloat_pos_inf(nfloat_ptr res, gr_ctx_t ctx)
{
    if (!(NFLOAT_CTX_FLAGS(ctx) & NFLOAT_ALLOW_INF))
        return GR_DOMAIN;

    NFLOAT_EXP(res) = NFLOAT_EXP_POS_INF;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}

int
nfloat_neg_inf(nfloat_ptr res, gr_ctx_t ctx)
{
    if (!(NFLOAT_CTX_FLAGS(ctx) & NFLOAT_ALLOW_INF))
        return GR_DOMAIN;

    NFLOAT_EXP(res) = NFLOAT_EXP_NEG_INF;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}

int
nfloat_nan(nfloat_ptr res, gr_ctx_t ctx)
{
    if (!(NFLOAT_CTX_FLAGS(ctx) & NFLOAT_ALLOW_NAN))
        return GR_UNABLE;

    NFLOAT_EXP(res) = NFLOAT_EXP_NAN;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}

int
_nfloat_underflow(nfloat_ptr res, int sgnbit, gr_ctx_t ctx)
{
    if (!(NFLOAT_CTX_FLAGS(ctx) & NFLOAT_ALLOW_UNDERFLOW))
        return GR_UNABLE;

    return nfloat_zero(res, ctx);
}

int
_nfloat_overflow(nfloat_ptr res, int sgnbit, gr_ctx_t ctx)
{
    if (!(NFLOAT_CTX_FLAGS(ctx) & NFLOAT_ALLOW_INF))
        return GR_UNABLE;

    return sgnbit ? nfloat_neg_inf(res, ctx) : nfloat_pos_inf(res, ctx);
}

int
nfloat_set_ui(nfloat_ptr res, ulong x, gr_ctx_t ctx)
{
    if (x == 0)
    {
        return nfloat_zero(res, ctx);
    }
    else
    {
        slong n = NFLOAT_CTX_NLIMBS(ctx);
        slong norm = flint_clz(x);

        NFLOAT_EXP(res) = FLINT_BITS - norm;
        NFLOAT_SGNBIT(res) = 0;
        flint_mpn_zero(NFLOAT_D(res), n - 1);
        NFLOAT_D(res)[n - 1] = x << norm;
        return GR_SUCCESS;
    }
}

int
nfloat_set_si(nfloat_ptr res, slong x, gr_ctx_t ctx)
{
    if (x == 0)
    {
        return nfloat_zero(res, ctx);
    }
    else
    {
        ulong ux = (x > 0) ? x : -x;
        slong n = NFLOAT_CTX_NLIMBS(ctx);
        slong norm = flint_clz(ux);

        NFLOAT_EXP(res) = FLINT_BITS - norm;
        NFLOAT_SGNBIT(res) = (x < 0);
        flint_mpn_zero(NFLOAT_D(res), n - 1);
        NFLOAT_D(res)[n - 1] = ux << norm;
        return GR_SUCCESS;
    }
}

int
nfloat_set_fmpz(nfloat_ptr res, const fmpz_t x, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*x))
    {
        return nfloat_set_si(res, *x, ctx);
    }
    else
    {
        mpz_ptr z = COEFF_TO_PTR(*x);
        slong zn = z->_mp_size;

        if (zn > 0)
            return _nfloat_set_mpn_2exp(res, z->_mp_d, zn, zn * FLINT_BITS, 0, ctx);
        else
            return _nfloat_set_mpn_2exp(res, z->_mp_d, -zn, -zn * FLINT_BITS, 1, ctx);
    }
}

/* todo: fast code */
int
nfloat_set_fmpq(nfloat_ptr res, const fmpq_t v, gr_ctx_t ctx)
{
    arf_t t;
    int status;
    arf_init(t);
    arf_set_fmpq(t, v, NFLOAT_CTX_PREC(ctx), ARF_RND_DOWN);
    status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);
    return status;
}

/* todo: fast code */
int
nfloat_set_d(nfloat_ptr res, double x, gr_ctx_t ctx)
{
    arf_t t;
    int status;
    arf_init(t);
    arf_set_d(t, x);
    status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);
    return status;
}

int
nfloat_set_str(nfloat_ptr res, const char * x, gr_ctx_t ctx)
{
    int status;

    arb_t t;
    arb_init(t);

    if (!arb_set_str(t, x, NFLOAT_CTX_PREC(ctx) + 20))
    {
        arf_set_round(arb_midref(t), arb_midref(t), NFLOAT_CTX_PREC(ctx), ARF_RND_NEAR);
        status = nfloat_set_arf(res, arb_midref(t), ctx);
    }
    else
    {
        status = gr_generic_set_str_ring_exponents(res, x, ctx);
    }

    arb_clear(t);
    return status;
}

int
nfloat_set_other(nfloat_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_NFLOAT:
            {
                slong nlimbs, x_nlimbs;

                if (NFLOAT_IS_SPECIAL(x))
                {
                    if (NFLOAT_IS_ZERO(x))
                        return nfloat_zero(res, ctx);
                    if (NFLOAT_IS_POS_INF(x))
                        return nfloat_pos_inf(res, ctx);
                    if (NFLOAT_IS_NEG_INF(x))
                        return nfloat_neg_inf(res, ctx);
                    return nfloat_nan(res, ctx);
                }

                nlimbs = NFLOAT_CTX_NLIMBS(ctx);
                x_nlimbs = NFLOAT_CTX_NLIMBS(x_ctx);

                NFLOAT_EXP(res) = NFLOAT_EXP(x);
                NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x);

                if (nlimbs <= x_nlimbs)
                {
                    flint_mpn_copyi(NFLOAT_D(res), NFLOAT_D(x) + x_nlimbs - nlimbs, nlimbs);
                }
                else
                {
                    flint_mpn_zero(NFLOAT_D(res), nlimbs - x_nlimbs);
                    flint_mpn_copyi(NFLOAT_D(res) + nlimbs - x_nlimbs, NFLOAT_D(x), x_nlimbs);
                }

                return GR_SUCCESS;
            }

        case GR_CTX_FMPZ:
            return nfloat_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return nfloat_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_FLOAT_ARF:
            return nfloat_set_arf(res, x, ctx);

        case GR_CTX_RR_ARB:
            return nfloat_set_arf(res, arb_midref((arb_srcptr) x), ctx);

        default:
            {
                int status;
                arf_t t;

                gr_ctx_t arf_ctx;
                arf_init(t);

                gr_ctx_init_real_float_arf(arf_ctx, NFLOAT_CTX_PREC(ctx));
                status = gr_set_other(t, x, x_ctx, arf_ctx);
                if (status == GR_SUCCESS)
                    status = nfloat_set_arf(res, t, ctx);

                arf_clear(t);
                gr_ctx_clear(arf_ctx);
                return status;
            }
    }
}

int
nfloat_set_arf(nfloat_ptr res, const arf_t x, gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            return nfloat_zero(res, ctx);
        else if (arf_is_pos_inf(x))
            return nfloat_pos_inf(res, ctx);
        else if (arf_is_neg_inf(x))
            return nfloat_neg_inf(res, ctx);
        else
            return nfloat_nan(res, ctx);
    }
    else
    {
        slong exp, xn;
        nn_srcptr xp;
        int sgnbit;

        ARF_GET_MPN_READONLY(xp, xn, x);

        exp = ARF_EXP(x);
        sgnbit = ARF_SGNBIT(x);

        if (COEFF_IS_MPZ(exp) || exp < NFLOAT_MIN_EXP || exp > NFLOAT_MAX_EXP)
        {
            if (fmpz_sgn(ARF_EXPREF(x)) < 0)
                return _nfloat_underflow(res, sgnbit, ctx);
            else
                return _nfloat_overflow(res, sgnbit, ctx);
        }

        _nfloat_set_mpn_2exp(res, xp, xn, exp, sgnbit, ctx);
    }

    return GR_SUCCESS;
}

int
nfloat_get_arf(arf_t res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            arf_zero(res);
        else if (NFLOAT_IS_POS_INF(x))
            arf_pos_inf(res);
        else if (NFLOAT_IS_NEG_INF(x))
            arf_neg_inf(res);
        else
            arf_nan(res);
    }
    else
    {
        if (!LIMB_MSB_IS_SET(NFLOAT_D(x)[NFLOAT_CTX_NLIMBS(ctx) - 1]))
        {
            flint_printf("bad nfloat!\n");
            flint_abort();
        }

        arf_set_mpn(res, NFLOAT_D(x), NFLOAT_CTX_NLIMBS(ctx), NFLOAT_SGNBIT(x));
        arf_mul_2exp_si(res, res, NFLOAT_EXP(x) - FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx));
    }

    return GR_SUCCESS;
}

int
nfloat_neg(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (res != x)
    {
        slong i, n = NFLOAT_CTX_DATA_NLIMBS(ctx);

        for (i = 0; i < n; i++)
            NFLOAT_DATA(res)[i] = NFLOAT_DATA(x)[i];
    }

    if (NFLOAT_IS_SPECIAL(res))
    {
        if (NFLOAT_IS_POS_INF(res))
            NFLOAT_EXP(res) = NFLOAT_EXP_NEG_INF;
        else if (NFLOAT_IS_NEG_INF(res))
            NFLOAT_EXP(res) = NFLOAT_EXP_POS_INF;
    }
    else
    {
        NFLOAT_SGNBIT(res) ^= 1;
    }

    return GR_SUCCESS;
}

int
nfloat_abs(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (res != x)
    {
        slong i, n = NFLOAT_CTX_DATA_NLIMBS(ctx);

        for (i = 0; i < n; i++)
            NFLOAT_DATA(res)[i] = NFLOAT_DATA(x)[i];
    }

    if (NFLOAT_IS_SPECIAL(res))
    {
        if (NFLOAT_IS_NEG_INF(res))
            NFLOAT_EXP(res) = NFLOAT_EXP_POS_INF;
    }
    else
    {
        NFLOAT_SGNBIT(res) = 0;
    }

    return GR_SUCCESS;
}

int
_nfloat_cmp(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    int xs, ys, mc;
    slong xe, ye;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        {
            if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
                return 0;
            if (NFLOAT_IS_POS_INF(x)) return NFLOAT_IS_POS_INF(y) ? 0 : 1;
            if (NFLOAT_IS_NEG_INF(x)) return NFLOAT_IS_NEG_INF(y) ? 0 : -1;
            if (NFLOAT_IS_POS_INF(y)) return -1;
            if (NFLOAT_IS_NEG_INF(y)) return 1;
        }

        if (NFLOAT_IS_ZERO(x) && NFLOAT_IS_ZERO(y)) return 0;
        if (NFLOAT_IS_ZERO(x)) return NFLOAT_SGNBIT(y) ? 1 : -1;
        if (NFLOAT_IS_ZERO(y)) return NFLOAT_SGNBIT(x) ? -1 : 1;
    }

    xs = NFLOAT_SGNBIT(x);
    ys = NFLOAT_SGNBIT(y);

    if (xs != ys)
        return xs ? -1 : 1;

    xe = NFLOAT_EXP(x);
    ye = NFLOAT_EXP(y);

    if (xe != ye)
        return ((xe < ye) ^ xs) ? -1 : 1;

    mc = mpn_cmp(NFLOAT_D(x), NFLOAT_D(y), NFLOAT_CTX_NLIMBS(ctx));

    if (mc != 0)
        return ((mc < 0) ^ xs) ? -1 : 1;

    return 0;
}

int
_nfloat_cmpabs(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    int mc;
    slong xe, ye;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        {
            if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
                return 0;
            if (NFLOAT_IS_INF(x)) return NFLOAT_IS_INF(y) ? 0 : 1;
            if (NFLOAT_IS_INF(y)) return -1;
        }

        if (NFLOAT_IS_ZERO(x) && NFLOAT_IS_ZERO(y)) return 0;
        if (NFLOAT_IS_ZERO(x)) return -1;
        if (NFLOAT_IS_ZERO(y)) return 1;
    }

    xe = NFLOAT_EXP(x);
    ye = NFLOAT_EXP(y);

    if (xe != ye)
        return (xe < ye) ? -1 : 1;

    mc = mpn_cmp(NFLOAT_D(x), NFLOAT_D(y), NFLOAT_CTX_NLIMBS(ctx));

    if (mc != 0)
        return (mc < 0) ? -1 : 1;

    return 0;
}

int
nfloat_cmp(int * res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
        return GR_UNABLE;

    *res = _nfloat_cmp(x, y, ctx);
    return GR_SUCCESS;
}

int
nfloat_cmpabs(int * res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
        return GR_UNABLE;

    *res = _nfloat_cmpabs(x, y, ctx);
    return GR_SUCCESS;
}

int
nfloat_sgn(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);

        if (NFLOAT_IS_POS_INF(x))
            return nfloat_one(res, ctx);

        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_neg_one(res, ctx);

        return nfloat_nan(res, ctx);
    }
    else
    {
        return NFLOAT_SGNBIT(x) ? nfloat_neg_one(res, ctx) : nfloat_one(res, ctx);
    }
}

int
nfloat_im(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    return nfloat_zero(res, ctx);
}

static mp_limb_pair_t
_flint_mpn_mulhigh_normalised2(nn_ptr rp, nn_srcptr xp, nn_srcptr yp, slong n)
{
    mp_limb_pair_t ret;

    FLINT_ASSERT(n >= 1);

    if (rp == xp || rp == yp)
    {
        nn_ptr t;
        TMP_INIT;
        TMP_START;
        t = TMP_ALLOC(sizeof(ulong) * n);
        ret = _flint_mpn_mulhigh_normalised2(t, xp, yp, n);
        flint_mpn_copyi(rp, t, n);
        TMP_END;
        return ret;
    }

    if (xp == yp)
        ret.m1 = flint_mpn_sqrhigh(rp, xp, n);
    else
        ret.m1 = flint_mpn_mulhigh_n(rp, xp, yp, n);

    if (LIMB_MSB_IS_SET(rp[n - 1]))
    {
        ret.m2 = 0;
    }
    else
    {
        ret.m2 = 1;
        mpn_lshift(rp, rp, n, 1);
        rp[0] |= (ret.m1 >> (FLINT_BITS - 1));
        ret.m1 <<= 1;
    }

    return ret;
}

/* handles aliasing */
FLINT_FORCE_INLINE
mp_limb_pair_t flint_mpn_mulhigh_normalised2(nn_ptr rp, nn_srcptr xp, nn_srcptr yp, slong n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_MULHIGH_NORMALISED_FUNC(n))
        return flint_mpn_mulhigh_normalised_func_tab[n](rp, xp, yp);
    else
        return _flint_mpn_mulhigh_normalised2(rp, xp, yp, n);
}

/* handles aliasing */
FLINT_FORCE_INLINE
mp_limb_pair_t flint_mpn_sqrhigh_normalised2(nn_ptr rp, nn_srcptr xp, slong n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_SQRHIGH_NORMALISED_FUNC(n))
        return flint_mpn_sqrhigh_normalised_func_tab[n](rp, xp);
    else
        return _flint_mpn_mulhigh_normalised2(rp, xp, xp, n);
}

FLINT_FORCE_INLINE int
n_signed_sub(ulong * r, ulong x, ulong y)
{
    if (x >= y)
    {
        *r = x - y;
        return 0;
    }
    else
    {
        *r = y - x;
        return 1;
    }
}

int
_nfloat_add_1(nfloat_ptr res, ulong x0, slong xexp, int xsgnbit, ulong y0, slong delta, gr_ctx_t ctx)
{
    ulong hi, lo;

    NFLOAT_SGNBIT(res) = xsgnbit;

    if (delta < FLINT_BITS)
    {
        add_ssaaaa(hi, lo, 0, x0, 0, y0 >> delta);

        if (hi == 0)
        {
            NFLOAT_D(res)[0] = lo;
            NFLOAT_EXP(res) = xexp;
        }
        else
        {
            NFLOAT_D(res)[0] = (lo >> 1) | (UWORD(1) << (FLINT_BITS - 1));
            NFLOAT_EXP(res) = xexp + 1;
            NFLOAT_HANDLE_OVERFLOW(res, ctx);
        }
    }
    else
    {
        NFLOAT_D(res)[0] = x0;
        NFLOAT_EXP(res) = xexp;
        NFLOAT_SGNBIT(res) = xsgnbit;
    }

    return GR_SUCCESS;
}

int
_nfloat_add_2(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
{
    ulong x1, x0, y1, y0, s2, s1, s0;

    NFLOAT_SGNBIT(res) = xsgnbit;

    x0 = x[0];
    x1 = x[1];

    if (delta < 2 * FLINT_BITS)
    {
        y0 = y[0];
        y1 = y[1];

        if (delta == 0)
            add_sssaaaaaa(s2, s1, s0, 0, x1, x0, 0, y1, y0);
        else if (delta < FLINT_BITS)
            add_sssaaaaaa(s2, s1, s0, 0, x1, x0, 0, y1 >> delta, (y0 >> delta) | (y1 << (FLINT_BITS - delta)));
        else
            add_sssaaaaaa(s2, s1, s0, 0, x1, x0, 0, 0, y1 >> (delta - FLINT_BITS));

        if (s2 == 0)
        {
            NFLOAT_D(res)[0] = s0;
            NFLOAT_D(res)[1] = s1;
            NFLOAT_EXP(res) = xexp;
        }
        else
        {
            NFLOAT_D(res)[0] = (s0 >> 1) | (s1 << (FLINT_BITS - 1));
            NFLOAT_D(res)[1] = (s1 >> 1) | (UWORD(1) << (FLINT_BITS - 1));
            NFLOAT_EXP(res) = xexp + 1;
            NFLOAT_HANDLE_OVERFLOW(res, ctx);
        }
    }
    else
    {
        NFLOAT_D(res)[0] = x0;
        NFLOAT_D(res)[1] = x1;
        NFLOAT_EXP(res) = xexp;
        NFLOAT_SGNBIT(res) = xsgnbit;
    }

    return GR_SUCCESS;
}

int
_nfloat_add_3(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
{
    ulong x2, x1, x0, y2, y1, y0, s3, s2, s1, s0;

    NFLOAT_SGNBIT(res) = xsgnbit;

    x0 = x[0];
    x1 = x[1];
    x2 = x[2];

    if (delta < FLINT_BITS)
    {
        y0 = y[0];
        y1 = y[1];
        y2 = y[2];

        if (delta != 0)
        {
            y0 = (y0 >> delta) | (y1 << (FLINT_BITS - delta));
            y1 = (y1 >> delta) | (y2 << (FLINT_BITS - delta));
            y2 = (y2 >> delta);
        }
    }
    else if (delta < 2 * FLINT_BITS)
    {
        delta -= FLINT_BITS;
        y0 = y[1];
        y1 = y[2];
        y2 = 0;

        if (delta != 0)
        {
            y0 = (y0 >> delta) | (y1 << (FLINT_BITS - delta));
            y1 = (y1 >> delta);
        }
    }
    else if (delta < 3 * FLINT_BITS)
    {
        delta -= 2 * FLINT_BITS;
        y0 = y[2] >> delta;
        y1 = 0;
        y2 = 0;
    }
    else
    {
        y0 = 0;
        y1 = 0;
        y2 = 0;
    }

    add_ssssaaaaaaaa(s3, s2, s1, s0, 0, x2, x1, x0, 0, y2, y1, y0);

    if (s3 == 0)
    {
        NFLOAT_D(res)[0] = s0;
        NFLOAT_D(res)[1] = s1;
        NFLOAT_D(res)[2] = s2;
        NFLOAT_EXP(res) = xexp;
    }
    else
    {
        NFLOAT_D(res)[0] = (s0 >> 1) | (s1 << (FLINT_BITS - 1));
        NFLOAT_D(res)[1] = (s1 >> 1) | (s2 << (FLINT_BITS - 1));
        NFLOAT_D(res)[2] = (s2 >> 1) | (UWORD(1) << (FLINT_BITS - 1));
        NFLOAT_EXP(res) = xexp + 1;
        NFLOAT_HANDLE_OVERFLOW(res, ctx);
    }

    return GR_SUCCESS;
}

int
_nfloat_add_4(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
{
    ulong x3, x2, x1, x0, y3, y2, y1, y0, s4, s3, s2, s1, s0;

    NFLOAT_SGNBIT(res) = xsgnbit;

    x0 = x[0];
    x1 = x[1];
    x2 = x[2];
    x3 = x[3];

    if (delta < FLINT_BITS)
    {
        y0 = y[0];
        y1 = y[1];
        y2 = y[2];
        y3 = y[3];

        if (delta != 0)
        {
            y0 = (y0 >> delta) | (y1 << (FLINT_BITS - delta));
            y1 = (y1 >> delta) | (y2 << (FLINT_BITS - delta));
            y2 = (y2 >> delta) | (y3 << (FLINT_BITS - delta));
            y3 = (y3 >> delta);
        }
    }
    else if (delta < 2 * FLINT_BITS)
    {
        delta -= FLINT_BITS;
        y0 = y[1];
        y1 = y[2];
        y2 = y[3];
        y3 = 0;

        if (delta != 0)
        {
            y0 = (y0 >> delta) | (y1 << (FLINT_BITS - delta));
            y1 = (y1 >> delta) | (y2 << (FLINT_BITS - delta));
            y2 = (y2 >> delta);
        }
    }
    else if (delta < 3 * FLINT_BITS)
    {
        delta -= 2 * FLINT_BITS;
        y0 = y[2];
        y1 = y[3];
        y2 = 0;
        y3 = 0;

        if (delta != 0)
        {
            y0 = (y0 >> delta) | (y1 << (FLINT_BITS - delta));
            y1 = (y1 >> delta);
        }
    }
    else if (delta < 4 * FLINT_BITS)
    {
        delta -= 3 * FLINT_BITS;
        y0 = y[3] >> delta;
        y1 = 0;
        y2 = 0;
        y3 = 0;
    }
    else
    {
        y0 = 0;
        y1 = 0;
        y2 = 0;
        y3 = 0;
    }

    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, 0, x3, x2, x1, x0, 0, y3, y2, y1, y0);

    if (s4 == 0)
    {
        NFLOAT_D(res)[0] = s0;
        NFLOAT_D(res)[1] = s1;
        NFLOAT_D(res)[2] = s2;
        NFLOAT_D(res)[3] = s3;
        NFLOAT_EXP(res) = xexp;
    }
    else
    {
        NFLOAT_D(res)[0] = (s0 >> 1) | (s1 << (FLINT_BITS - 1));
        NFLOAT_D(res)[1] = (s1 >> 1) | (s2 << (FLINT_BITS - 1));
        NFLOAT_D(res)[2] = (s2 >> 1) | (s3 << (FLINT_BITS - 1));
        NFLOAT_D(res)[3] = (s3 >> 1) | (UWORD(1) << (FLINT_BITS - 1));
        NFLOAT_EXP(res) = xexp + 1;
        NFLOAT_HANDLE_OVERFLOW(res, ctx);
    }

    return GR_SUCCESS;
}

int
_nfloat_sub_1(nfloat_ptr res, ulong x0, slong xexp, int xsgnbit, ulong y0, slong delta, gr_ctx_t ctx)
{
    ulong u, u0;
    slong norm = 0;

    if (delta <= 1)
    {
        if (delta == 0)
        {
            xsgnbit ^= n_signed_sub(&u, x0, y0);

            if (u == 0)
                return nfloat_zero(res, ctx);

            norm = flint_clz(u);
            u <<= norm;
            xexp -= norm;
        }
        else
        {
            sub_ddmmss(u, u0, x0, 0, y0 >> delta, y0 << (FLINT_BITS - delta));

            if (FLINT_UNLIKELY(u == 0))
            {
                u = u0;
                xexp -= FLINT_BITS;
            }
            else if (!LIMB_MSB_IS_SET(u))
            {
                norm = flint_clz(u);
                u = (u << norm) | (u0 >> (FLINT_BITS - norm));
                xexp -= norm;
            }
        }
    }
    else if (delta < FLINT_BITS)
    {
        u = x0 - (y0 >> delta);
        NFLOAT_SGNBIT(res) = xsgnbit;

        if (!LIMB_MSB_IS_SET(u))
        {
            u <<= 1;
            xexp--;
        }
    }
    else
    {
        u = x0;
    }

    NFLOAT_EXP(res) = xexp;
    NFLOAT_SGNBIT(res) = xsgnbit;
    NFLOAT_D(res)[0] = u;
    NFLOAT_HANDLE_UNDERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
_nfloat_sub_2(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
{
    ulong x1, x0, y1, y0, s1, s0, sb;
    slong norm;

    NFLOAT_SGNBIT(res) = xsgnbit;

    x0 = x[0];
    x1 = x[1];

    y0 = y[0];
    y1 = y[1];

    if (delta <= 1)
    {
        if (delta == 0)
        {
            if (x1 > y1 || (x1 == y1 && x0 >= y0))
            {
                sub_ddmmss(s1, s0, x1, x0, y1, y0);

                if (s1 == 0 && s0 == 0)
                    return nfloat_zero(res, ctx);
            }
            else
            {
                sub_ddmmss(s1, s0, y1, y0, x1, x0);
                xsgnbit = !xsgnbit;
            }

            sb = 0;
        }
        else
        {
            sub_dddmmmsss(s1, s0, sb,
                          x1, x0, 0,
                          y1 >> 1,
                          (y0 >> 1) | (y1 << (FLINT_BITS - 1)),
                          y0 << (FLINT_BITS - 1));
        }

        if (FLINT_UNLIKELY(s1 == 0))
        {
            if (s0 == 0)
            {
                s1 = sb;
                s0 = 0;
                xexp -= 2 * FLINT_BITS;
            }
            else
            {
                s1 = s0;
                s0 = sb;
                sb = 0;
                xexp -= FLINT_BITS;
            }
        }

        if (!LIMB_MSB_IS_SET(s1))
        {
            norm = flint_clz(s1);
            s1 = (s1 << norm) | (s0 >> (FLINT_BITS - norm));
            s0 = (s0 << norm) | (sb >> (FLINT_BITS - norm));
            xexp -= norm;
        }
    }
    else
    {
        if (delta < FLINT_BITS)
        {
            sub_ddmmss(s1, s0, x1, x0, y1 >> delta, (y0 >> delta) | (y1 << (FLINT_BITS - delta)));
        }
        else if (delta < 2 * FLINT_BITS)
        {
            delta -= FLINT_BITS;
            sub_ddmmss(s1, s0, x1, x0, 0, y1 >> delta);
        }
        else
        {
            s0 = x0;
            s1 = x1;
        }

        if (!LIMB_MSB_IS_SET(s1))
        {
            s1 = (s1 << 1) | (s0 >> (FLINT_BITS - 1));
            s0 = (s0 << 1);
            xexp--;
        }
    }

    NFLOAT_EXP(res) = xexp;
    NFLOAT_SGNBIT(res) = xsgnbit;
    NFLOAT_D(res)[0] = s0;
    NFLOAT_D(res)[1] = s1;
    NFLOAT_HANDLE_UNDERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
_nfloat_sub_3(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
{
    ulong x2, x1, x0, y2, y1, y0, s2, s1, s0, sb;
    slong norm;

    NFLOAT_SGNBIT(res) = xsgnbit;

    x0 = x[0];
    x1 = x[1];
    x2 = x[2];

    if (delta <= 1)
    {
        y0 = y[0];
        y1 = y[1];
        y2 = y[2];

        if (delta == 0)
        {
            if (x2 > y2 || (x2 == y2 && (x1 > y1 || (x1 == y1 && x0 >= y0))))
            {
                sub_dddmmmsss(s2, s1, s0, x2, x1, x0, y2, y1, y0);

                if (s2 == 0 && s1 == 0 && s0 == 0)
                    return nfloat_zero(res, ctx);
            }
            else
            {
                sub_dddmmmsss(s2, s1, s0, y2, y1, y0, x2, x1, x0);
                xsgnbit = !xsgnbit;
            }

            sb = 0;
        }
        else
        {
            sub_ddddmmmmssss(s2, s1, s0, sb,
                             x2, x1, x0, 0,
                             y2 >> 1,
                             (y1 >> 1) | (y2 << (FLINT_BITS - 1)),
                             (y0 >> 1) | (y1 << (FLINT_BITS - 1)),
                             y0 << (FLINT_BITS - 1));
        }

        if (FLINT_UNLIKELY(s2 == 0))
        {
            if (s1 == 0)
            {
                if (s0 == 0)
                {
                    s2 = sb;
                    s1 = 0;
                    s0 = 0;
                    sb = 0;
                    xexp -= 3 * FLINT_BITS;
                }
                else
                {
                    s2 = s0;
                    s1 = sb;
                    s0 = 0;
                    sb = 0;
                    xexp -= 2 * FLINT_BITS;
                }
            }
            else
            {
                s2 = s1;
                s1 = s0;
                s0 = sb;
                sb = 0;
                xexp -= FLINT_BITS;
            }
        }

        if (!LIMB_MSB_IS_SET(s2))
        {
            norm = flint_clz(s2);
            s2 = (s2 << norm) | (s1 >> (FLINT_BITS - norm));
            s1 = (s1 << norm) | (s0 >> (FLINT_BITS - norm));
            s0 = (s0 << norm) | (sb >> (FLINT_BITS - norm));
            xexp -= norm;
        }
    }
    else
    {
        if (delta < FLINT_BITS)
        {
            y0 = y[0];
            y1 = y[1];
            y2 = y[2];

            y0 = (y0 >> delta) | (y1 << (FLINT_BITS - delta));
            y1 = (y1 >> delta) | (y2 << (FLINT_BITS - delta));
            y2 = (y2 >> delta);
        }
        else if (delta < 2 * FLINT_BITS)
        {
            delta -= FLINT_BITS;
            y0 = y[1];
            y1 = y[2];
            y2 = 0;

            if (delta != 0)
            {
                y0 = (y0 >> delta) | (y1 << (FLINT_BITS - delta));
                y1 = (y1 >> delta);
            }
        }
        else if (delta < 3 * FLINT_BITS)
        {
            delta -= 2 * FLINT_BITS;
            y0 = y[2] >> delta;
            y1 = 0;
            y2 = 0;
        }
        else
        {
            y0 = 0;
            y1 = 0;
            y2 = 0;
        }

        sub_dddmmmsss(s2, s1, s0, x2, x1, x0, y2, y1, y0);

        if (!LIMB_MSB_IS_SET(s2))
        {
            s2 = (s2 << 1) | (s1 >> (FLINT_BITS - 1));
            s1 = (s1 << 1) | (s0 >> (FLINT_BITS - 1));
            s0 = (s0 << 1);
            xexp--;
        }
    }

    NFLOAT_EXP(res) = xexp;
    NFLOAT_SGNBIT(res) = xsgnbit;
    NFLOAT_D(res)[0] = s0;
    NFLOAT_D(res)[1] = s1;
    NFLOAT_D(res)[2] = s2;
    NFLOAT_HANDLE_UNDERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
_nfloat_sub_4(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
{
    ulong x3, x2, x1, x0, y3, y2, y1, y0, s3, s2, s1, s0, sb;
    slong norm;

    NFLOAT_SGNBIT(res) = xsgnbit;

    x0 = x[0];
    x1 = x[1];
    x2 = x[2];
    x3 = x[3];

    if (delta <= 1)
    {
        y0 = y[0];
        y1 = y[1];
        y2 = y[2];
        y3 = y[3];

        if (delta == 0)
        {
            if (x3 > y3 || (x3 == y3 && (x2 > y2 || (x2 == y2 && (x1 > y1 || (x1 == y1 && x0 >= y0))))))
            {
                sub_ddddmmmmssss(s3, s2, s1, s0, x3, x2, x1, x0, y3, y2, y1, y0);

                if (s3 == 0 && s2 == 0 && s1 == 0 && s0 == 0)
                    return nfloat_zero(res, ctx);
            }
            else
            {
                sub_ddddmmmmssss(s3, s2, s1, s0, y3, y2, y1, y0, x3, x2, x1, x0);
                xsgnbit = !xsgnbit;
            }

            sb = 0;
        }
        else
        {
            sub_dddddmmmmmsssss(s3, s2, s1, s0, sb,
                                x3, x2, x1, x0, 0,
                                y3 >> 1,
                                (y2 >> 1) | (y3 << (FLINT_BITS - 1)),
                                (y1 >> 1) | (y2 << (FLINT_BITS - 1)),
                                (y0 >> 1) | (y1 << (FLINT_BITS - 1)),
                                y0 << (FLINT_BITS - 1));
        }

        if (FLINT_UNLIKELY(s3 == 0))
        {
            if (s2 == 0)
            {
                if (s1 == 0)
                {
                    if (s0 == 0)
                    {
                        s3 = sb;
                        s2 = 0;
                        s1 = 0;
                        s0 = 0;
                        sb = 0;
                        xexp -= 4 * FLINT_BITS;
                    }
                    else
                    {
                        s3 = s0;
                        s2 = sb;
                        s1 = 0;
                        s0 = 0;
                        sb = 0;
                        xexp -= 3 * FLINT_BITS;
                    }
                }
                else
                {
                    s3 = s1;
                    s2 = s0;
                    s1 = sb;
                    s0 = 0;
                    sb = 0;
                    xexp -= 2 * FLINT_BITS;
                }
            }
            else
            {
                s3 = s2;
                s2 = s1;
                s1 = s0;
                s0 = sb;
                sb = 0;
                xexp -= FLINT_BITS;
            }
        }

        if (!LIMB_MSB_IS_SET(s3))
        {
            norm = flint_clz(s3);
            s3 = (s3 << norm) | (s2 >> (FLINT_BITS - norm));
            s2 = (s2 << norm) | (s1 >> (FLINT_BITS - norm));
            s1 = (s1 << norm) | (s0 >> (FLINT_BITS - norm));
            s0 = (s0 << norm) | (sb >> (FLINT_BITS - norm));
            xexp -= norm;
        }
    }
    else
    {
        if (delta < FLINT_BITS)
        {
            y0 = y[0];
            y1 = y[1];
            y2 = y[2];
            y3 = y[3];

            y0 = (y0 >> delta) | (y1 << (FLINT_BITS - delta));
            y1 = (y1 >> delta) | (y2 << (FLINT_BITS - delta));
            y2 = (y2 >> delta) | (y3 << (FLINT_BITS - delta));
            y3 = (y3 >> delta);
        }
        else if (delta < 2 * FLINT_BITS)
        {
            delta -= FLINT_BITS;
            y0 = y[1];
            y1 = y[2];
            y2 = y[3];
            y3 = 0;

            if (delta != 0)
            {
                y0 = (y0 >> delta) | (y1 << (FLINT_BITS - delta));
                y1 = (y1 >> delta) | (y2 << (FLINT_BITS - delta));
                y2 = (y2 >> delta);
            }
        }
        else if (delta < 3 * FLINT_BITS)
        {
            delta -= 2 * FLINT_BITS;
            y0 = y[2];
            y1 = y[3];
            y2 = 0;
            y3 = 0;

            if (delta != 0)
            {
                y0 = (y0 >> delta) | (y1 << (FLINT_BITS - delta));
                y1 = (y1 >> delta);
            }
        }
        else if (delta < 4 * FLINT_BITS)
        {
            delta -= 3 * FLINT_BITS;
            y0 = y[3] >> delta;
            y1 = 0;
            y2 = 0;
            y3 = 0;
        }
        else
        {
            y0 = 0;
            y1 = 0;
            y2 = 0;
            y3 = 0;
        }

        sub_ddddmmmmssss(s3, s2, s1, s0, x3, x2, x1, x0, y3, y2, y1, y0);

        if (!LIMB_MSB_IS_SET(s3))
        {
            s3 = (s3 << 1) | (s2 >> (FLINT_BITS - 1));
            s2 = (s2 << 1) | (s1 >> (FLINT_BITS - 1));
            s1 = (s1 << 1) | (s0 >> (FLINT_BITS - 1));
            s0 = (s0 << 1);
            xexp--;
        }
    }

    NFLOAT_EXP(res) = xexp;
    NFLOAT_SGNBIT(res) = xsgnbit;
    NFLOAT_D(res)[0] = s0;
    NFLOAT_D(res)[1] = s1;
    NFLOAT_D(res)[2] = s2;
    NFLOAT_D(res)[3] = s3;
    NFLOAT_HANDLE_UNDERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
_nfloat_add_n(nfloat_ptr res, nn_srcptr xd, slong xexp, int xsgnbit, nn_srcptr yd, slong delta, slong nlimbs, gr_ctx_t ctx)
{
    slong shift_limbs, shift_bits;
    ulong cy;
    ulong t[NFLOAT_MAX_LIMBS];

    NFLOAT_SGNBIT(res) = xsgnbit;

    if (delta == 0)
    {
        cy = mpn_add_n(NFLOAT_D(res), xd, yd, nlimbs);
    }
    else if (delta < FLINT_BITS)
    {
        mpn_rshift(t, yd, nlimbs, delta);
        cy = mpn_add_n(NFLOAT_D(res), xd, t, nlimbs);
    }
    else if (delta < nlimbs * FLINT_BITS)
    {
        shift_limbs = delta / FLINT_BITS;
        shift_bits = delta % FLINT_BITS;
        if (shift_bits == 0)
            flint_mpn_copyi(t, yd + shift_limbs, nlimbs - shift_limbs);
        else
            mpn_rshift(t, yd + shift_limbs, nlimbs - shift_limbs, shift_bits);
        cy = mpn_add(NFLOAT_D(res), xd, nlimbs, t, nlimbs - shift_limbs);
    }
    else
    {
        flint_mpn_copyi(NFLOAT_D(res), xd, nlimbs);
        cy = 0;
    }

    if (cy == 0)
    {
        NFLOAT_EXP(res) = xexp;
    }
    else
    {
        mpn_rshift(NFLOAT_D(res), NFLOAT_D(res), nlimbs, 1);
        NFLOAT_D(res)[nlimbs - 1] |= (UWORD(1) << (FLINT_BITS - 1));
        NFLOAT_EXP(res) = xexp + 1;
        NFLOAT_HANDLE_OVERFLOW(res, ctx);
    }

    return GR_SUCCESS;
}

int
_nfloat_sub_n(nfloat_ptr res, nn_srcptr xd, slong xexp, int xsgnbit, nn_srcptr yd, slong delta, slong nlimbs, gr_ctx_t ctx)
{
    slong shift_limbs, shift_bits, n, norm;
    ulong t[NFLOAT_MAX_LIMBS];
    ulong sb;

    /* In some experiments,

            ~10% of additions have shift == 0
            ~10% of additions have shift == 1
            ~60% of additions have 2 <= shift < FLINT_BITS
            ~20% of additions have shift >= FLINT_BITS

        The shift == 1 case is the only case other than delta == 0 where
        we can have significant cancellation. In this case, we need to
        subtract the shifted-out bit of y to guarantee that the output
        approximates the exact difference with high relative accuracy.
    */
    if (delta <= 1)
    {
        if (delta == 0)
        {
            xsgnbit ^= flint_mpn_signed_sub_n(NFLOAT_D(res), xd, yd, nlimbs);
            sb = 0;
        }
        else
        {
            mpn_rshift(t, yd, nlimbs, 1);
            sb = yd[0] << (FLINT_BITS - 1);
            mpn_sub_n(NFLOAT_D(res), xd, t, nlimbs);
            mpn_sub_1(NFLOAT_D(res), NFLOAT_D(res), nlimbs, sb != 0);
        }

        n = nlimbs;
        MPN_NORM(NFLOAT_D(res), n);

        /* We have at least one full limb of cancellation */
        if (FLINT_UNLIKELY(n != nlimbs))
        {
            if (n == 0)
            {
                if (sb == 0)
                {
                    return nfloat_zero(res, ctx);
                }
                else
                {
                    NFLOAT_D(res)[nlimbs - 1] = sb;
                    flint_mpn_zero(NFLOAT_D(res), nlimbs - 1);
                    xexp -= nlimbs * FLINT_BITS;
                }
            }
            else
            {
                /* The copy is avoidable when we also have a shift,
                   but since this is an unlikely, branch, it's
                   not really worth optimizing. */
                flint_mpn_copyd(NFLOAT_D(res) + nlimbs - n, NFLOAT_D(res), n);
                flint_mpn_zero(NFLOAT_D(res), nlimbs - n - 1);
                NFLOAT_D(res)[nlimbs - n - 1] = sb;
                sb = 0;
                xexp -= (nlimbs - n) * FLINT_BITS;
            }
        }

        norm = flint_clz(NFLOAT_D(res)[nlimbs - 1]);

        if (norm)
        {
            mpn_lshift(NFLOAT_D(res), NFLOAT_D(res), nlimbs, norm);
            NFLOAT_D(res)[0] |= (sb >> (FLINT_BITS - norm));
            xexp -= norm;
        }
    }
    else
    {
        if (delta < FLINT_BITS)
        {
            mpn_rshift(t, yd, nlimbs, delta);
            mpn_sub_n(NFLOAT_D(res), xd, t, nlimbs);
        }
        else if (delta < nlimbs * FLINT_BITS)
        {
            shift_limbs = delta / FLINT_BITS;
            shift_bits = delta % FLINT_BITS;
            if (shift_bits == 0)
                flint_mpn_copyi(t, yd + shift_limbs, nlimbs - shift_limbs);
            else
                mpn_rshift(t, yd + shift_limbs, nlimbs - shift_limbs, shift_bits);
            mpn_sub(NFLOAT_D(res), xd, nlimbs, t, nlimbs - shift_limbs);
        }
        else
        {
            flint_mpn_copyi(NFLOAT_D(res), xd, nlimbs);
        }

        if (!LIMB_MSB_IS_SET(NFLOAT_D(res)[nlimbs - 1]))
        {
            mpn_lshift(NFLOAT_D(res), NFLOAT_D(res), nlimbs, 1);
            xexp--;
        }
    }

    NFLOAT_SGNBIT(res) = xsgnbit;
    NFLOAT_EXP(res) = xexp;
    NFLOAT_HANDLE_UNDERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
nfloat_add(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong xexp, yexp, delta, nlimbs;
    int xsgnbit, ysgnbit;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_set(res, y, ctx);
        if (NFLOAT_IS_ZERO(y))
            return nfloat_set(res, x, ctx);
        if (NFLOAT_IS_POS_INF(x) && NFLOAT_IS_POS_INF(y))
            return nfloat_pos_inf(res, ctx);
        if (NFLOAT_IS_NEG_INF(x) && NFLOAT_IS_NEG_INF(y))
            return nfloat_neg_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    xexp = NFLOAT_EXP(x);
    yexp = NFLOAT_EXP(y);

    if (xexp < yexp)
    {
        FLINT_SWAP(nfloat_srcptr, x, y);
        FLINT_SWAP(slong, xexp, yexp);
    }

    xsgnbit = NFLOAT_SGNBIT(x);
    ysgnbit = NFLOAT_SGNBIT(y);

    delta = xexp - yexp;

    if (xsgnbit == ysgnbit)
    {
        if (nlimbs == 1)
            return _nfloat_add_1(res, NFLOAT_D(x)[0], xexp, xsgnbit, NFLOAT_D(y)[0], delta, ctx);
        else if (nlimbs == 2)
            return _nfloat_add_2(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else if (nlimbs == 3)
            return _nfloat_add_3(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else if (nlimbs == 4)
            return _nfloat_add_4(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else
            return _nfloat_add_n(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, nlimbs, ctx);
    }
    else
    {
        if (nlimbs == 1)
            return _nfloat_sub_1(res, NFLOAT_D(x)[0], xexp, xsgnbit, NFLOAT_D(y)[0], delta, ctx);
        else if (nlimbs == 2)
            return _nfloat_sub_2(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else if (nlimbs == 3)
            return _nfloat_sub_3(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else if (nlimbs == 4)
            return _nfloat_sub_4(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else
            return _nfloat_sub_n(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, nlimbs, ctx);
    }
}

int
nfloat_sub(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong xexp, yexp, delta, nlimbs;
    int xsgnbit, ysgnbit;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_neg(res, y, ctx);
        if (NFLOAT_IS_ZERO(y))
            return nfloat_set(res, x, ctx);
        if (NFLOAT_IS_POS_INF(x) && NFLOAT_IS_NEG_INF(y))
            return nfloat_pos_inf(res, ctx);
        if (NFLOAT_IS_NEG_INF(x) && NFLOAT_IS_POS_INF(y))
            return nfloat_neg_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    xexp = NFLOAT_EXP(x);
    yexp = NFLOAT_EXP(y);

    xsgnbit = NFLOAT_SGNBIT(x);
    ysgnbit = NFLOAT_SGNBIT(y) ^ 1;

    if (xexp < yexp)
    {
        FLINT_SWAP(nfloat_srcptr, x, y);
        FLINT_SWAP(slong, xexp, yexp);
        FLINT_SWAP(int, xsgnbit, ysgnbit);
    }

    delta = xexp - yexp;

    if (xsgnbit == ysgnbit)
    {
        if (nlimbs == 1)
            return _nfloat_add_1(res, NFLOAT_D(x)[0], xexp, xsgnbit, NFLOAT_D(y)[0], delta, ctx);
        else if (nlimbs == 2)
            return _nfloat_add_2(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else if (nlimbs == 3)
            return _nfloat_add_3(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else if (nlimbs == 4)
            return _nfloat_add_4(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else
            return _nfloat_add_n(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, nlimbs, ctx);
    }
    else
    {
        if (nlimbs == 1)
            return _nfloat_sub_1(res, NFLOAT_D(x)[0], xexp, xsgnbit, NFLOAT_D(y)[0], delta, ctx);
        else if (nlimbs == 2)
            return _nfloat_sub_2(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else if (nlimbs == 3)
            return _nfloat_sub_3(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else if (nlimbs == 4)
            return _nfloat_sub_4(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, ctx);
        else
            return _nfloat_sub_n(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, nlimbs, ctx);
    }
}

int
_nfloat_vec_aors_1(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, int subtract, slong len, gr_ctx_t ctx)
{
    slong i, stride;
    slong xexp, yexp, delta;
    int xsgnbit, ysgnbit;
    nfloat_srcptr xi, yi;
    nfloat_ptr ri;
    int status = GR_SUCCESS;

    stride = NFLOAT_HEADER_LIMBS + 1;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        yi = (nn_srcptr) y + i * stride;
        ri = (nn_ptr) res + i * stride;

        xexp = NFLOAT_EXP(xi);
        yexp = NFLOAT_EXP(yi);

        if (yexp == NFLOAT_EXP_ZERO)
        {
            flint_mpn_copyi(ri, xi, NFLOAT_HEADER_LIMBS + 1);
            continue;
        }

        if (xexp == NFLOAT_EXP_ZERO)
        {
            flint_mpn_copyi(ri, yi, NFLOAT_HEADER_LIMBS + 1);
            if (subtract)
                NFLOAT_SGNBIT(ri) = !NFLOAT_SGNBIT(ri);
            continue;
        }

        xsgnbit = NFLOAT_SGNBIT(xi);
        ysgnbit = NFLOAT_SGNBIT(yi) ^ subtract;

        delta = xexp - yexp;

        if (xsgnbit == ysgnbit)
        {
            if (delta >= 0)
                status |= _nfloat_add_1(ri, NFLOAT_D(xi)[0], xexp, xsgnbit, NFLOAT_D(yi)[0], delta, ctx);
            else
                status |= _nfloat_add_1(ri, NFLOAT_D(yi)[0], yexp, ysgnbit, NFLOAT_D(xi)[0], -delta, ctx);
        }
        else
        {
            if (delta >= 0)
                status |= _nfloat_sub_1(ri, NFLOAT_D(xi)[0], xexp, xsgnbit, NFLOAT_D(yi)[0], delta, ctx);
            else
                status |= _nfloat_sub_1(ri, NFLOAT_D(yi)[0], yexp, ysgnbit, NFLOAT_D(xi)[0], -delta, ctx);
        }
    }

    return status;
}

int
_nfloat_vec_aors_2(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, int subtract, slong len, gr_ctx_t ctx)
{
    slong i, stride;
    slong xexp, yexp, delta;
    int xsgnbit, ysgnbit;
    nfloat_srcptr xi, yi;
    nfloat_ptr ri;
    int status = GR_SUCCESS;

    stride = NFLOAT_HEADER_LIMBS + 2;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        yi = (nn_srcptr) y + i * stride;
        ri = (nn_ptr) res + i * stride;

        xexp = NFLOAT_EXP(xi);
        yexp = NFLOAT_EXP(yi);

        if (yexp == NFLOAT_EXP_ZERO)
        {
            flint_mpn_copyi(ri, xi, NFLOAT_HEADER_LIMBS + 2);
            continue;
        }

        if (xexp == NFLOAT_EXP_ZERO)
        {
            flint_mpn_copyi(ri, yi, NFLOAT_HEADER_LIMBS + 2);
            if (subtract)
                NFLOAT_SGNBIT(ri) = !NFLOAT_SGNBIT(ri);
            continue;
        }

        xsgnbit = NFLOAT_SGNBIT(xi);
        ysgnbit = NFLOAT_SGNBIT(yi) ^ subtract;

        delta = xexp - yexp;

        if (xsgnbit == ysgnbit)
        {
            if (delta >= 0)
                status |= _nfloat_add_2(ri, NFLOAT_D(xi), xexp, xsgnbit, NFLOAT_D(yi), delta, ctx);
            else
                status |= _nfloat_add_2(ri, NFLOAT_D(yi), yexp, ysgnbit, NFLOAT_D(xi), -delta, ctx);
        }
        else
        {
            if (delta >= 0)
                status |= _nfloat_sub_2(ri, NFLOAT_D(xi), xexp, xsgnbit, NFLOAT_D(yi), delta, ctx);
            else
                status |= _nfloat_sub_2(ri, NFLOAT_D(yi), yexp, ysgnbit, NFLOAT_D(xi), -delta, ctx);
        }
    }

    return status;
}

int
_nfloat_vec_aors_3(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, int subtract, slong len, gr_ctx_t ctx)
{
    slong i, stride;
    slong xexp, yexp, delta;
    int xsgnbit, ysgnbit;
    nfloat_srcptr xi, yi;
    nfloat_ptr ri;
    int status = GR_SUCCESS;

    stride = NFLOAT_HEADER_LIMBS + 3;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        yi = (nn_srcptr) y + i * stride;
        ri = (nn_ptr) res + i * stride;

        xexp = NFLOAT_EXP(xi);
        yexp = NFLOAT_EXP(yi);

        if (yexp == NFLOAT_EXP_ZERO)
        {
            flint_mpn_copyi(ri, xi, NFLOAT_HEADER_LIMBS + 3);
            continue;
        }

        if (xexp == NFLOAT_EXP_ZERO)
        {
            flint_mpn_copyi(ri, yi, NFLOAT_HEADER_LIMBS + 3);
            if (subtract)
                NFLOAT_SGNBIT(ri) = !NFLOAT_SGNBIT(ri);
            continue;
        }

        xsgnbit = NFLOAT_SGNBIT(xi);
        ysgnbit = NFLOAT_SGNBIT(yi) ^ subtract;

        delta = xexp - yexp;

        if (xsgnbit == ysgnbit)
        {
            if (delta >= 0)
                status |= _nfloat_add_3(ri, NFLOAT_D(xi), xexp, xsgnbit, NFLOAT_D(yi), delta, ctx);
            else
                status |= _nfloat_add_3(ri, NFLOAT_D(yi), yexp, ysgnbit, NFLOAT_D(xi), -delta, ctx);
        }
        else
        {
            if (delta >= 0)
                status |= _nfloat_sub_3(ri, NFLOAT_D(xi), xexp, xsgnbit, NFLOAT_D(yi), delta, ctx);
            else
                status |= _nfloat_sub_3(ri, NFLOAT_D(yi), yexp, ysgnbit, NFLOAT_D(xi), -delta, ctx);
        }
    }

    return status;
}

int
_nfloat_vec_aors_4(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, int subtract, slong len, gr_ctx_t ctx)
{
    slong i, stride;
    slong xexp, yexp, delta;
    int xsgnbit, ysgnbit;
    nfloat_srcptr xi, yi;
    nfloat_ptr ri;
    int status = GR_SUCCESS;

    stride = NFLOAT_HEADER_LIMBS + 4;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        yi = (nn_srcptr) y + i * stride;
        ri = (nn_ptr) res + i * stride;

        xexp = NFLOAT_EXP(xi);
        yexp = NFLOAT_EXP(yi);

        if (yexp == NFLOAT_EXP_ZERO)
        {
            flint_mpn_copyi(ri, xi, NFLOAT_HEADER_LIMBS + 4);
            continue;
        }

        if (xexp == NFLOAT_EXP_ZERO)
        {
            flint_mpn_copyi(ri, yi, NFLOAT_HEADER_LIMBS + 4);
            if (subtract)
                NFLOAT_SGNBIT(ri) = !NFLOAT_SGNBIT(ri);
            continue;
        }

        xsgnbit = NFLOAT_SGNBIT(xi);
        ysgnbit = NFLOAT_SGNBIT(yi) ^ subtract;

        delta = xexp - yexp;

        if (xsgnbit == ysgnbit)
        {
            if (delta >= 0)
                status |= _nfloat_add_4(ri, NFLOAT_D(xi), xexp, xsgnbit, NFLOAT_D(yi), delta, ctx);
            else
                status |= _nfloat_add_4(ri, NFLOAT_D(yi), yexp, ysgnbit, NFLOAT_D(xi), -delta, ctx);
        }
        else
        {
            if (delta >= 0)
                status |= _nfloat_sub_4(ri, NFLOAT_D(xi), xexp, xsgnbit, NFLOAT_D(yi), delta, ctx);
            else
                status |= _nfloat_sub_4(ri, NFLOAT_D(yi), yexp, ysgnbit, NFLOAT_D(xi), -delta, ctx);
        }
    }

    return status;
}

int
_nfloat_vec_add(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
{
    slong i, stride, nlimbs;
    nfloat_srcptr xi, yi;
    nfloat_ptr ri;
    int status = GR_SUCCESS;

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    if (!NFLOAT_CTX_HAS_INF_NAN(ctx))
    {
        if (nlimbs == 1)
            return _nfloat_vec_aors_1(res, x, y, 0, len, ctx);
        if (nlimbs == 2)
            return _nfloat_vec_aors_2(res, x, y, 0, len, ctx);
        if (nlimbs == 3)
            return _nfloat_vec_aors_3(res, x, y, 0, len, ctx);
        if (nlimbs == 4)
            return _nfloat_vec_aors_4(res, x, y, 0, len, ctx);
    }

    stride = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        yi = (nn_srcptr) y + i * stride;
        ri = (nn_ptr) res + i * stride;

        status |= nfloat_add(ri, xi, yi, ctx);
    }

    return status;
}

int
_nfloat_vec_sub(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
{
    slong i, stride, nlimbs;
    nfloat_srcptr xi, yi;
    nfloat_ptr ri;
    int status = GR_SUCCESS;

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    if (!NFLOAT_CTX_HAS_INF_NAN(ctx))
    {
        if (nlimbs == 1)
            return _nfloat_vec_aors_1(res, x, y, 1, len, ctx);
        if (nlimbs == 2)
            return _nfloat_vec_aors_2(res, x, y, 1, len, ctx);
        if (nlimbs == 3)
            return _nfloat_vec_aors_3(res, x, y, 1, len, ctx);
        if (nlimbs == 4)
            return _nfloat_vec_aors_4(res, x, y, 1, len, ctx);
    }

    stride = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        yi = (nn_srcptr) y + i * stride;
        ri = (nn_ptr) res + i * stride;

        status |= nfloat_sub(ri, xi, yi, ctx);
    }

    return status;
}

int
nfloat_mul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    mp_limb_pair_t mul_res;
    slong nlimbs;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_CTX_FLAGS(ctx) & (NFLOAT_ALLOW_INF | NFLOAT_ALLOW_NAN))
        {
            if (NFLOAT_IS_ZERO(x) && !NFLOAT_IS_SPECIAL(y))
                return nfloat_zero(res, ctx);
            if (NFLOAT_IS_ZERO(y) && !NFLOAT_IS_SPECIAL(x))
                return nfloat_zero(res, ctx);
            return GR_UNABLE;
        }
        else
        {
            return nfloat_zero(res, ctx);
        }
    }

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    if (nlimbs == 1)
    {
        ulong hi, lo;

        umul_ppmm(hi, lo, NFLOAT_D(x)[0], NFLOAT_D(y)[0]);

        if (LIMB_MSB_IS_SET(hi))
        {
            NFLOAT_D(res)[0] = hi;
            NFLOAT_EXP(res) = NFLOAT_EXP(x) + NFLOAT_EXP(y);
        }
        else
        {
            NFLOAT_D(res)[0] = (hi << 1) | (lo >> (FLINT_BITS - 1));
            NFLOAT_EXP(res) = NFLOAT_EXP(x) + NFLOAT_EXP(y) - 1;
        }
    }
    else if (nlimbs == 2)
    {
        ulong r3, r2, r1, FLINT_SET_BUT_UNUSED(r0);

        FLINT_MPN_MUL_2X2(r3, r2, r1, r0, NFLOAT_D(x)[1], NFLOAT_D(x)[0], NFLOAT_D(y)[1], NFLOAT_D(y)[0]);

        if (LIMB_MSB_IS_SET(r3))
        {
            NFLOAT_D(res)[0] = r2;
            NFLOAT_D(res)[1] = r3;
            NFLOAT_EXP(res) = NFLOAT_EXP(x) + NFLOAT_EXP(y);
        }
        else
        {
            NFLOAT_D(res)[0] = (r2 << 1) | (r1 >> (FLINT_BITS - 1));
            NFLOAT_D(res)[1] = (r3 << 1) | (r2 >> (FLINT_BITS - 1));
            NFLOAT_EXP(res) = NFLOAT_EXP(x) + NFLOAT_EXP(y) - 1;
        }
    }
    else
    {
        mul_res = flint_mpn_mulhigh_normalised2(NFLOAT_D(res), NFLOAT_D(x), NFLOAT_D(y), NFLOAT_CTX_NLIMBS(ctx));
        NFLOAT_EXP(res) = NFLOAT_EXP(x) + NFLOAT_EXP(y) - mul_res.m2;
    }

    NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x) ^ NFLOAT_SGNBIT(y);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
nfloat_sqr(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    mp_limb_pair_t mul_res;
    slong nlimbs;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        else
            return nfloat_abs(res, x, ctx);
    }

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    if (nlimbs == 1)
    {
        ulong hi, lo;

        umul_ppmm(hi, lo, NFLOAT_D(x)[0], NFLOAT_D(x)[0]);

        if (LIMB_MSB_IS_SET(hi))
        {
            NFLOAT_D(res)[0] = hi;
            NFLOAT_EXP(res) = 2 * NFLOAT_EXP(x);
        }
        else
        {
            NFLOAT_D(res)[0] = (hi << 1) | (lo >> (FLINT_BITS - 1));
            NFLOAT_EXP(res) = 2 * NFLOAT_EXP(x) - 1;
        }
    }
    else if (nlimbs == 2)
    {
        ulong r3, r2, r1, FLINT_SET_BUT_UNUSED(r0);

        FLINT_MPN_SQR_2X2(r3, r2, r1, r0, NFLOAT_D(x)[1], NFLOAT_D(x)[0]);

        if (LIMB_MSB_IS_SET(r3))
        {
            NFLOAT_D(res)[0] = r2;
            NFLOAT_D(res)[1] = r3;
            NFLOAT_EXP(res) = 2 * NFLOAT_EXP(x);
        }
        else
        {
            NFLOAT_D(res)[0] = (r2 << 1) | (r1 >> (FLINT_BITS - 1));
            NFLOAT_D(res)[1] = (r3 << 1) | (r2 >> (FLINT_BITS - 1));
            NFLOAT_EXP(res) = 2 * NFLOAT_EXP(x) - 1;
        }
    }
    else
    {
        mul_res = flint_mpn_sqrhigh_normalised2(NFLOAT_D(res), NFLOAT_D(x), NFLOAT_CTX_NLIMBS(ctx));
        NFLOAT_EXP(res) = 2 * NFLOAT_EXP(x) - mul_res.m2;
    }

    NFLOAT_SGNBIT(res) = 0;
    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
_nfloat_vec_mul_1(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
{
    slong i, stride;
    slong xexp, yexp;
    nfloat_srcptr xi, yi;
    nfloat_ptr ri;
    ulong hi, lo;
    int xsgnbit, ysgnbit;
    int status = GR_SUCCESS;

    stride = NFLOAT_HEADER_LIMBS + 1;

    if (x == y)
    {
        for (i = 0; i < len; i++)
        {
            xi = (nn_srcptr) x + i * stride;
            ri = (nn_ptr) res + i * stride;

            xexp = NFLOAT_EXP(xi);

            if (xexp == NFLOAT_EXP_ZERO)
            {
                NFLOAT_EXP(ri) = NFLOAT_EXP_ZERO;
                NFLOAT_SGNBIT(ri) = 0;
                continue;
            }

            umul_ppmm(hi, lo, NFLOAT_D(xi)[0], NFLOAT_D(xi)[0]);

            if (LIMB_MSB_IS_SET(hi))
            {
                NFLOAT_D(ri)[0] = hi;
                NFLOAT_EXP(ri) = 2 * xexp;
            }
            else
            {
                NFLOAT_D(ri)[0] = (hi << 1) | (lo >> (FLINT_BITS - 1));
                NFLOAT_EXP(ri) = 2 * xexp - 1;
            }

            NFLOAT_SGNBIT(ri) = 0;

            if (FLINT_UNLIKELY(NFLOAT_EXP(ri) < NFLOAT_MIN_EXP))
                status |= _nfloat_underflow(ri, NFLOAT_SGNBIT(ri), ctx);
            if (FLINT_UNLIKELY(NFLOAT_EXP(ri) > NFLOAT_MAX_EXP))
                status |= _nfloat_overflow(ri, NFLOAT_SGNBIT(ri), ctx);
        }
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            xi = (nn_srcptr) x + i * stride;
            yi = (nn_srcptr) y + i * stride;
            ri = (nn_ptr) res + i * stride;

            xexp = NFLOAT_EXP(xi);
            yexp = NFLOAT_EXP(yi);

            if (xexp == NFLOAT_EXP_ZERO || yexp == NFLOAT_EXP_ZERO)
            {
                NFLOAT_EXP(ri) = NFLOAT_EXP_ZERO;
                NFLOAT_SGNBIT(ri) = 0;
                continue;
            }

            umul_ppmm(hi, lo, NFLOAT_D(xi)[0], NFLOAT_D(yi)[0]);

            xsgnbit = NFLOAT_SGNBIT(xi);
            ysgnbit = NFLOAT_SGNBIT(yi);

            if (LIMB_MSB_IS_SET(hi))
            {
                NFLOAT_D(ri)[0] = hi;
                NFLOAT_EXP(ri) = xexp + yexp;
            }
            else
            {
                NFLOAT_D(ri)[0] = (hi << 1) | (lo >> (FLINT_BITS - 1));
                NFLOAT_EXP(ri) = xexp + yexp - 1;
            }

            NFLOAT_SGNBIT(ri) = xsgnbit ^ ysgnbit;

            if (FLINT_UNLIKELY(NFLOAT_EXP(ri) < NFLOAT_MIN_EXP))
                status |= _nfloat_underflow(ri, NFLOAT_SGNBIT(ri), ctx);
            if (FLINT_UNLIKELY(NFLOAT_EXP(ri) > NFLOAT_MAX_EXP))
                status |= _nfloat_overflow(ri, NFLOAT_SGNBIT(ri), ctx);
        }
    }

    return status;
}

int
_nfloat_vec_mul_2(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
{
    slong i, stride;
    slong xexp, yexp;
    nfloat_srcptr xi, yi;
    nfloat_ptr ri;
    ulong r3, r2, r1, FLINT_SET_BUT_UNUSED(r0);
    int xsgnbit, ysgnbit;
    int status = GR_SUCCESS;

    stride = NFLOAT_HEADER_LIMBS + 2;

    if (x == y)
    {
        for (i = 0; i < len; i++)
        {
            xi = (nn_srcptr) x + i * stride;
            ri = (nn_ptr) res + i * stride;

            xexp = NFLOAT_EXP(xi);

            if (xexp == NFLOAT_EXP_ZERO)
            {
                NFLOAT_EXP(ri) = NFLOAT_EXP_ZERO;
                NFLOAT_SGNBIT(ri) = 0;
                continue;
            }

            FLINT_MPN_SQR_2X2(r3, r2, r1, r0, NFLOAT_D(xi)[1], NFLOAT_D(xi)[0]);

            if (LIMB_MSB_IS_SET(r3))
            {
                NFLOAT_D(ri)[0] = r2;
                NFLOAT_D(ri)[1] = r3;
                NFLOAT_EXP(ri) = 2 * xexp;
            }
            else
            {
                NFLOAT_D(ri)[0] = (r2 << 1) | (r1 >> (FLINT_BITS - 1));
                NFLOAT_D(ri)[1] = (r3 << 1) | (r2 >> (FLINT_BITS - 1));
                NFLOAT_EXP(ri) = 2 * xexp - 1;
            }

            NFLOAT_SGNBIT(ri) = 0;

            if (FLINT_UNLIKELY(NFLOAT_EXP(ri) < NFLOAT_MIN_EXP))
                status |= _nfloat_underflow(ri, NFLOAT_SGNBIT(ri), ctx);
            if (FLINT_UNLIKELY(NFLOAT_EXP(ri) > NFLOAT_MAX_EXP))
                status |= _nfloat_overflow(ri, NFLOAT_SGNBIT(ri), ctx);
        }
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            xi = (nn_srcptr) x + i * stride;
            yi = (nn_srcptr) y + i * stride;
            ri = (nn_ptr) res + i * stride;

            xexp = NFLOAT_EXP(xi);
            yexp = NFLOAT_EXP(yi);

            if (xexp == NFLOAT_EXP_ZERO || yexp == NFLOAT_EXP_ZERO)
            {
                NFLOAT_EXP(ri) = NFLOAT_EXP_ZERO;
                NFLOAT_SGNBIT(ri) = 0;
                continue;
            }

            FLINT_MPN_MUL_2X2(r3, r2, r1, r0, NFLOAT_D(xi)[1], NFLOAT_D(xi)[0], NFLOAT_D(yi)[1], NFLOAT_D(yi)[0]);

            xsgnbit = NFLOAT_SGNBIT(xi);
            ysgnbit = NFLOAT_SGNBIT(yi);

            if (LIMB_MSB_IS_SET(r3))
            {
                NFLOAT_D(ri)[0] = r2;
                NFLOAT_D(ri)[1] = r3;
                NFLOAT_EXP(ri) = xexp + yexp;
            }
            else
            {
                NFLOAT_D(ri)[0] = (r2 << 1) | (r1 >> (FLINT_BITS - 1));
                NFLOAT_D(ri)[1] = (r3 << 1) | (r2 >> (FLINT_BITS - 1));
                NFLOAT_EXP(ri) = xexp + yexp - 1;
            }

            NFLOAT_SGNBIT(ri) = xsgnbit ^ ysgnbit;

            if (FLINT_UNLIKELY(NFLOAT_EXP(ri) < NFLOAT_MIN_EXP))
                status |= _nfloat_underflow(ri, NFLOAT_SGNBIT(ri), ctx);
            if (FLINT_UNLIKELY(NFLOAT_EXP(ri) > NFLOAT_MAX_EXP))
                status |= _nfloat_overflow(ri, NFLOAT_SGNBIT(ri), ctx);
        }
    }

    return status;
}

int
_nfloat_vec_mul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
{
    slong i, stride, nlimbs;
    nfloat_srcptr xi, yi;
    nfloat_ptr ri;
    int status = GR_SUCCESS;

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    if (!NFLOAT_CTX_HAS_INF_NAN(ctx))
    {
        if (nlimbs == 1)
            return _nfloat_vec_mul_1(res, x, y, len, ctx);
        if (nlimbs == 2)
            return _nfloat_vec_mul_2(res, x, y, len, ctx);
    }

    stride = NFLOAT_CTX_DATA_NLIMBS(ctx);

    if (x == y)
    {
        for (i = 0; i < len; i++)
        {
            xi = (nn_srcptr) x + i * stride;
            ri = (nn_ptr) res + i * stride;

            status |= nfloat_sqr(ri, xi, ctx);
        }
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            xi = (nn_srcptr) x + i * stride;
            yi = (nn_srcptr) y + i * stride;
            ri = (nn_ptr) res + i * stride;

            status |= nfloat_mul(ri, xi, yi, ctx);
        }
    }

    return status;
}

int
_nfloat_vec_mul_scalar_1(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong i, stride;
    slong xexp, yexp;
    nfloat_srcptr xi;
    nfloat_ptr ri;
    ulong y0, hi, lo;
    int xsgnbit, ysgnbit;
    int status = GR_SUCCESS;

    yexp = NFLOAT_EXP(y);
    y0 = NFLOAT_D(y)[0];
    ysgnbit = NFLOAT_SGNBIT(y);

    if (yexp == NFLOAT_EXP_ZERO)
        return _nfloat_vec_zero(res, len, ctx);

    stride = NFLOAT_HEADER_LIMBS + 1;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        ri = (nn_ptr) res + i * stride;

        xexp = NFLOAT_EXP(xi);

        if (xexp == NFLOAT_EXP_ZERO)
        {
            NFLOAT_EXP(ri) = NFLOAT_EXP_ZERO;
            NFLOAT_SGNBIT(ri) = 0;
            continue;
        }

        umul_ppmm(hi, lo, NFLOAT_D(xi)[0], y0);

        xsgnbit = NFLOAT_SGNBIT(xi);

        if (LIMB_MSB_IS_SET(hi))
        {
            NFLOAT_D(ri)[0] = hi;
            NFLOAT_EXP(ri) = xexp + yexp;
        }
        else
        {
            NFLOAT_D(ri)[0] = (hi << 1) | (lo >> (FLINT_BITS - 1));
            NFLOAT_EXP(ri) = xexp + yexp - 1;
        }

        NFLOAT_SGNBIT(ri) = xsgnbit ^ ysgnbit;

        if (FLINT_UNLIKELY(NFLOAT_EXP(ri) < NFLOAT_MIN_EXP))
            status |= _nfloat_underflow(ri, NFLOAT_SGNBIT(ri), ctx);
        if (FLINT_UNLIKELY(NFLOAT_EXP(ri) > NFLOAT_MAX_EXP))
            status |= _nfloat_overflow(ri, NFLOAT_SGNBIT(ri), ctx);
    }

    return status;
}

int
_nfloat_vec_mul_scalar_2(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong i, stride;
    slong xexp, yexp;
    nfloat_srcptr xi;
    nfloat_ptr ri;
    ulong y0, y1;
    ulong r3, r2, r1, FLINT_SET_BUT_UNUSED(r0);
    int xsgnbit, ysgnbit;
    int status = GR_SUCCESS;

    yexp = NFLOAT_EXP(y);
    y0 = NFLOAT_D(y)[0];
    y1 = NFLOAT_D(y)[1];
    ysgnbit = NFLOAT_SGNBIT(y);

    if (yexp == NFLOAT_EXP_ZERO)
        return _nfloat_vec_zero(res, len, ctx);

    stride = NFLOAT_HEADER_LIMBS + 2;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        ri = (nn_ptr) res + i * stride;

        xexp = NFLOAT_EXP(xi);

        if (xexp == NFLOAT_EXP_ZERO)
        {
            NFLOAT_EXP(ri) = NFLOAT_EXP_ZERO;
            NFLOAT_SGNBIT(ri) = 0;
            continue;
        }

        FLINT_MPN_MUL_2X2(r3, r2, r1, r0, NFLOAT_D(xi)[1], NFLOAT_D(xi)[0], y1, y0);

        xsgnbit = NFLOAT_SGNBIT(xi);

        if (LIMB_MSB_IS_SET(r3))
        {
            NFLOAT_D(ri)[0] = r2;
            NFLOAT_D(ri)[1] = r3;
            NFLOAT_EXP(ri) = xexp + yexp;
        }
        else
        {
            NFLOAT_D(ri)[0] = (r2 << 1) | (r1 >> (FLINT_BITS - 1));
            NFLOAT_D(ri)[1] = (r3 << 1) | (r2 >> (FLINT_BITS - 1));
            NFLOAT_EXP(ri) = xexp + yexp - 1;
        }

        NFLOAT_SGNBIT(ri) = xsgnbit ^ ysgnbit;

        if (FLINT_UNLIKELY(NFLOAT_EXP(ri) < NFLOAT_MIN_EXP))
            status |= _nfloat_underflow(ri, NFLOAT_SGNBIT(ri), ctx);
        if (FLINT_UNLIKELY(NFLOAT_EXP(ri) > NFLOAT_MAX_EXP))
            status |= _nfloat_overflow(ri, NFLOAT_SGNBIT(ri), ctx);
    }

    return status;
}

int
_nfloat_vec_mul_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong i, stride, nlimbs;
    nfloat_srcptr xi;
    nfloat_ptr ri;
    int status = GR_SUCCESS;

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    if (!NFLOAT_CTX_HAS_INF_NAN(ctx))
    {
        if (nlimbs == 1)
            return _nfloat_vec_mul_scalar_1(res, x, len, y, ctx);
        if (nlimbs == 2)
            return _nfloat_vec_mul_scalar_2(res, x, len, y, ctx);
    }

    stride = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        ri = (nn_ptr) res + i * stride;

        status |= nfloat_mul(ri, xi, y, ctx);
    }

    return status;
}
int
nfloat_mul_2exp_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        return nfloat_set(res, x, ctx);
    }
    else
    {
        /* todo */
        if (y < NFLOAT_MIN_EXP || y > NFLOAT_MAX_EXP)
            return GR_UNABLE;

        nfloat_set(res, x, ctx);
        NFLOAT_EXP(res) += y;
        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
        return GR_SUCCESS;
    }
}

int
nfloat_addmul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    ulong t[NFLOAT_MAX_ALLOC];
    int status;

    status = nfloat_mul(t, x, y, ctx);
    status |= nfloat_add(res, res, t, ctx);
    return status;
}

int
nfloat_submul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    ulong t[NFLOAT_MAX_ALLOC];
    int status;

    status = nfloat_mul(t, x, y, ctx);
    status |= nfloat_sub(res, res, t, ctx);
    return status;
}

/* No infs or nans allowed by context; y != 0 */
int
_nfloat_vec_aorsmul_scalar_1(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, int subtract, gr_ctx_t ctx)
{
    slong yexp, rexp, texp;
    slong i, stride;
    nfloat_srcptr xi;
    nfloat_ptr ri;
    ulong t[NFLOAT_HEADER_LIMBS + 1];
    ulong y0, hi, lo;
    int rsgnbit, tsgnbit, ysgnbit;
    int status = GR_SUCCESS;

    yexp = NFLOAT_EXP(y);
    ysgnbit = NFLOAT_SGNBIT(y) ^ subtract;
    y0 = NFLOAT_D(y)[0];

    stride = NFLOAT_HEADER_LIMBS + 1;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        ri = (nn_ptr) res + i * stride;

        if (NFLOAT_IS_ZERO(xi))
            continue;

        /* status |= nfloat_mul(t, xi, y, ctx); */
        umul_ppmm(hi, lo, NFLOAT_D(xi)[0], y0);

        if (LIMB_MSB_IS_SET(hi))
        {
            NFLOAT_D(t)[0] = hi;
            NFLOAT_EXP(t) = NFLOAT_EXP(xi) + yexp;
        }
        else
        {
            NFLOAT_D(t)[0] = (hi << 1) | (lo >> (FLINT_BITS - 1));
            NFLOAT_EXP(t) = NFLOAT_EXP(xi) + yexp - 1;
        }

        NFLOAT_SGNBIT(t) = NFLOAT_SGNBIT(xi) ^ ysgnbit;
        texp = NFLOAT_EXP(t);

        /* The product could have underflowed or overflowed. By assumption that
           we have no infs or nans, quit. */
        if (FLINT_UNLIKELY(texp < NFLOAT_MIN_EXP || texp > NFLOAT_MAX_EXP))
            return GR_UNABLE;

        if (NFLOAT_IS_ZERO(ri))
        {
            flint_mpn_copyi(ri, t, NFLOAT_HEADER_LIMBS + 1);
            continue;
        }

        rexp = NFLOAT_EXP(ri);
        tsgnbit = NFLOAT_SGNBIT(t);
        rsgnbit = NFLOAT_SGNBIT(ri);

        if (rexp >= texp)
            status |= (NFLOAT_SGNBIT(t) == NFLOAT_SGNBIT(ri) ? _nfloat_add_1 : _nfloat_sub_1)(ri, NFLOAT_D(ri)[0], rexp, rsgnbit, NFLOAT_D(t)[0], rexp - texp, ctx);
        else
            status |= (NFLOAT_SGNBIT(t) == NFLOAT_SGNBIT(ri) ? _nfloat_add_1 : _nfloat_sub_1)(ri, NFLOAT_D(t)[0], texp, tsgnbit, NFLOAT_D(ri)[0], texp - rexp, ctx);
    }

    return status;
}

/* No infs or nans allowed by context; y != 0 */
int
_nfloat_vec_aorsmul_scalar_2(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, int subtract, gr_ctx_t ctx)
{
    slong yexp, rexp, texp;
    slong i, stride;
    nfloat_srcptr xi;
    nfloat_ptr ri;
    ulong t[NFLOAT_HEADER_LIMBS + 2];
    ulong y0, y1, r3, r2, r1, FLINT_SET_BUT_UNUSED(r0);
    int rsgnbit, tsgnbit, ysgnbit;
    int status = GR_SUCCESS;

    yexp = NFLOAT_EXP(y);
    ysgnbit = NFLOAT_SGNBIT(y) ^ subtract;
    y0 = NFLOAT_D(y)[0];
    y1 = NFLOAT_D(y)[1];

    stride = NFLOAT_HEADER_LIMBS + 2;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        ri = (nn_ptr) res + i * stride;

        if (NFLOAT_IS_ZERO(xi))
            continue;

        /* status |= nfloat_mul(t, xi, y, ctx); */

        FLINT_MPN_MUL_2X2(r3, r2, r1, r0, NFLOAT_D(xi)[1], NFLOAT_D(xi)[0], y1, y0);

        if (LIMB_MSB_IS_SET(r3))
        {
            NFLOAT_D(t)[0] = r2;
            NFLOAT_D(t)[1] = r3;
            NFLOAT_EXP(t) = NFLOAT_EXP(xi) + yexp;
        }
        else
        {
            NFLOAT_D(t)[0] = (r2 << 1) | (r1 >> (FLINT_BITS - 1));
            NFLOAT_D(t)[1] = (r3 << 1) | (r2 >> (FLINT_BITS - 1));
            NFLOAT_EXP(t) = NFLOAT_EXP(xi) + yexp - 1;
        }

        NFLOAT_SGNBIT(t) = NFLOAT_SGNBIT(xi) ^ ysgnbit;
        texp = NFLOAT_EXP(t);

        /* The product could have underflowed or overflowed. By assumption that
           we have no infs or nans, quit. */
        if (FLINT_UNLIKELY(texp < NFLOAT_MIN_EXP || texp > NFLOAT_MAX_EXP))
            return GR_UNABLE;

        if (NFLOAT_IS_ZERO(ri))
        {
            flint_mpn_copyi(ri, t, NFLOAT_HEADER_LIMBS + 2);
            continue;
        }

        rexp = NFLOAT_EXP(ri);
        tsgnbit = NFLOAT_SGNBIT(t);
        rsgnbit = NFLOAT_SGNBIT(ri);

        if (rexp >= texp)
            status |= (NFLOAT_SGNBIT(t) == NFLOAT_SGNBIT(ri) ? _nfloat_add_2 : _nfloat_sub_2)(ri, NFLOAT_D(ri), rexp, rsgnbit, NFLOAT_D(t), rexp - texp, ctx);
        else
            status |= (NFLOAT_SGNBIT(t) == NFLOAT_SGNBIT(ri) ? _nfloat_add_2 : _nfloat_sub_2)(ri, NFLOAT_D(t), texp, tsgnbit, NFLOAT_D(ri), texp - rexp, ctx);
    }

    return status;
}

/* No infs or nans allowed by context; y != 0 */
int
_nfloat_vec_aorsmul_scalar_3(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, int subtract, gr_ctx_t ctx)
{
    slong yexp, rexp, texp;
    slong i, stride;
    nfloat_srcptr xi;
    nfloat_ptr ri;
    ulong t[NFLOAT_HEADER_LIMBS + 3];
    int rsgnbit, tsgnbit, ysgnbit;
    int status = GR_SUCCESS;
    mp_limb_pair_t mul_res;

    yexp = NFLOAT_EXP(y);
    ysgnbit = NFLOAT_SGNBIT(y) ^ subtract;

    stride = NFLOAT_HEADER_LIMBS + 3;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        ri = (nn_ptr) res + i * stride;

        if (NFLOAT_IS_ZERO(xi))
            continue;

        /* status |= nfloat_mul(t, xi, y, ctx); */
        mul_res = flint_mpn_mulhigh_normalised2(NFLOAT_D(t), NFLOAT_D(xi), NFLOAT_D(y), 3);
        NFLOAT_EXP(t) = NFLOAT_EXP(xi) + yexp - mul_res.m2;
        NFLOAT_SGNBIT(t) = NFLOAT_SGNBIT(xi) ^ ysgnbit;
        texp = NFLOAT_EXP(t);

        /* The product could have underflowed or overflowed. By assumption that
           we have no infs or nans, quit. */
        if (FLINT_UNLIKELY(texp < NFLOAT_MIN_EXP || texp > NFLOAT_MAX_EXP))
            return GR_UNABLE;

        if (NFLOAT_IS_ZERO(ri))
        {
            flint_mpn_copyi(ri, t, NFLOAT_HEADER_LIMBS + 3);
            continue;
        }

        rexp = NFLOAT_EXP(ri);
        tsgnbit = NFLOAT_SGNBIT(t);
        rsgnbit = NFLOAT_SGNBIT(ri);

        if (rexp >= texp)
            status |= (NFLOAT_SGNBIT(t) == NFLOAT_SGNBIT(ri) ? _nfloat_add_3 : _nfloat_sub_3)(ri, NFLOAT_D(ri), rexp, rsgnbit, NFLOAT_D(t), rexp - texp, ctx);
        else
            status |= (NFLOAT_SGNBIT(t) == NFLOAT_SGNBIT(ri) ? _nfloat_add_3 : _nfloat_sub_3)(ri, NFLOAT_D(t), texp, tsgnbit, NFLOAT_D(ri), texp - rexp, ctx);
    }

    return status;
}

int
_nfloat_vec_aorsmul_scalar_4(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, int subtract, gr_ctx_t ctx)
{
    slong yexp, rexp, texp;
    slong i, stride;
    nfloat_srcptr xi;
    nfloat_ptr ri;
    ulong t[NFLOAT_HEADER_LIMBS + 4];
    int rsgnbit, tsgnbit, ysgnbit;
    int status = GR_SUCCESS;
    mp_limb_pair_t mul_res;

    yexp = NFLOAT_EXP(y);
    ysgnbit = NFLOAT_SGNBIT(y) ^ subtract;

    stride = NFLOAT_HEADER_LIMBS + 4;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        ri = (nn_ptr) res + i * stride;

        if (NFLOAT_IS_ZERO(xi))
            continue;

        /* status |= nfloat_mul(t, xi, y, ctx); */
        mul_res = flint_mpn_mulhigh_normalised2(NFLOAT_D(t), NFLOAT_D(xi), NFLOAT_D(y), 4);
        NFLOAT_EXP(t) = NFLOAT_EXP(xi) + yexp - mul_res.m2;
        NFLOAT_SGNBIT(t) = NFLOAT_SGNBIT(xi) ^ ysgnbit;
        texp = NFLOAT_EXP(t);

        /* The product could have underflowed or overflowed. By assumption that
           we have no infs or nans, quit. */
        if (FLINT_UNLIKELY(texp < NFLOAT_MIN_EXP || texp > NFLOAT_MAX_EXP))
            return GR_UNABLE;

        if (NFLOAT_IS_ZERO(ri))
        {
            flint_mpn_copyi(ri, t, NFLOAT_HEADER_LIMBS + 4);
            continue;
        }

        rexp = NFLOAT_EXP(ri);
        tsgnbit = NFLOAT_SGNBIT(t);
        rsgnbit = NFLOAT_SGNBIT(ri);

        if (rexp >= texp)
            status |= (NFLOAT_SGNBIT(t) == NFLOAT_SGNBIT(ri) ? _nfloat_add_4 : _nfloat_sub_4)(ri, NFLOAT_D(ri), rexp, rsgnbit, NFLOAT_D(t), rexp - texp, ctx);
        else
            status |= (NFLOAT_SGNBIT(t) == NFLOAT_SGNBIT(ri) ? _nfloat_add_4 : _nfloat_sub_4)(ri, NFLOAT_D(t), texp, tsgnbit, NFLOAT_D(ri), texp - rexp, ctx);
    }

    return status;
}

int
_nfloat_vec_aorsmul_scalar_n(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, int subtract, slong nlimbs, gr_ctx_t ctx)
{
    slong yexp, rexp, texp;
    slong i, stride;
    nfloat_srcptr xi;
    nfloat_ptr ri;
    ulong t[NFLOAT_MAX_ALLOC];
    int rsgnbit, tsgnbit, ysgnbit;
    int status = GR_SUCCESS;
    mp_limb_pair_t mul_res;

    yexp = NFLOAT_EXP(y);
    ysgnbit = NFLOAT_SGNBIT(y) ^ subtract;

    stride = NFLOAT_HEADER_LIMBS + nlimbs;

    for (i = 0; i < len; i++)
    {
        xi = (nn_srcptr) x + i * stride;
        ri = (nn_ptr) res + i * stride;

        if (NFLOAT_IS_ZERO(xi))
            continue;

        /* status |= nfloat_mul(t, xi, y, ctx); */
        mul_res = flint_mpn_mulhigh_normalised2(NFLOAT_D(t), NFLOAT_D(xi), NFLOAT_D(y), nlimbs);
        NFLOAT_EXP(t) = NFLOAT_EXP(xi) + yexp - mul_res.m2;
        NFLOAT_SGNBIT(t) = NFLOAT_SGNBIT(xi) ^ ysgnbit;
        texp = NFLOAT_EXP(t);

        /* The product could have underflowed or overflowed. By assumption that
           we have no infs or nans, quit. */
        if (FLINT_UNLIKELY(texp < NFLOAT_MIN_EXP || texp > NFLOAT_MAX_EXP))
            return GR_UNABLE;

        if (NFLOAT_IS_ZERO(ri))
        {
            flint_mpn_copyi(ri, t, NFLOAT_HEADER_LIMBS + nlimbs);
            continue;
        }

        rexp = NFLOAT_EXP(ri);
        tsgnbit = NFLOAT_SGNBIT(t);
        rsgnbit = NFLOAT_SGNBIT(ri);

        if (rexp >= texp)
            status |= (NFLOAT_SGNBIT(t) == NFLOAT_SGNBIT(ri) ? _nfloat_add_n : _nfloat_sub_n)(ri, NFLOAT_D(ri), rexp, rsgnbit, NFLOAT_D(t), rexp - texp, nlimbs, ctx);
        else
            status |= (NFLOAT_SGNBIT(t) == NFLOAT_SGNBIT(ri) ? _nfloat_add_n : _nfloat_sub_n)(ri, NFLOAT_D(t), texp, tsgnbit, NFLOAT_D(ri), texp - rexp, nlimbs, ctx);
    }

    return status;
}

int
_nfloat_vec_addmul_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)
{
    if (!NFLOAT_CTX_HAS_INF_NAN(ctx))
    {
        slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);

        if (NFLOAT_IS_ZERO(y))
            return GR_SUCCESS;

        if (nlimbs == 1)
            return _nfloat_vec_aorsmul_scalar_1(res, x, len, y, 0, ctx);
        if (nlimbs == 2)
            return _nfloat_vec_aorsmul_scalar_2(res, x, len, y, 0, ctx);
        if (nlimbs == 3)
            return _nfloat_vec_aorsmul_scalar_3(res, x, len, y, 0, ctx);
        if (nlimbs == 4)
            return _nfloat_vec_aorsmul_scalar_4(res, x, len, y, 0, ctx);
        return _nfloat_vec_aorsmul_scalar_n(res, x, len, y, 0, nlimbs, ctx);
    }
    else
    {
        return gr_generic_vec_scalar_addmul(res, x, len, y, ctx);
    }
}

int
_nfloat_vec_submul_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)
{
    if (!NFLOAT_CTX_HAS_INF_NAN(ctx))
    {
        slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);

        if (NFLOAT_IS_ZERO(y))
            return GR_SUCCESS;

        if (nlimbs == 1)
            return _nfloat_vec_aorsmul_scalar_1(res, x, len, y, 1, ctx);
        if (nlimbs == 2)
            return _nfloat_vec_aorsmul_scalar_2(res, x, len, y, 1, ctx);
        if (nlimbs == 3)
            return _nfloat_vec_aorsmul_scalar_3(res, x, len, y, 1, ctx);
        if (nlimbs == 4)
            return _nfloat_vec_aorsmul_scalar_4(res, x, len, y, 1, ctx);
        return _nfloat_vec_aorsmul_scalar_n(res, x, len, y, 1, nlimbs, ctx);
    }
    else
    {
        return gr_generic_vec_scalar_submul(res, x, len, y, ctx);
    }
}

int
nfloat_inv(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    mpfr_t rf, xf;
    slong prec = NFLOAT_CTX_PREC(ctx);
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_nan(res, ctx);
        else
            return nfloat_nan(res, ctx); /* todo */
    }

    if (NFLOAT_D(x)[nlimbs - 1] == (UWORD(1) << (FLINT_BITS - 1)) &&
        flint_mpn_zero_p(NFLOAT_D(x), nlimbs - 1))
    {
        nfloat_set(res, x, ctx);
        NFLOAT_EXP(res) = 2 - NFLOAT_EXP(res);
        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
        return GR_SUCCESS;
    }

    /* todo: make sure aliasing is correct */

    rf->_mpfr_d = NFLOAT_D(res);
    rf->_mpfr_prec = prec;
    rf->_mpfr_sign = 1;
    rf->_mpfr_exp = 0;

    xf->_mpfr_d = NFLOAT_D(x);
    xf->_mpfr_prec = prec;
    xf->_mpfr_sign = 1;
    xf->_mpfr_exp = 0;

    mpfr_ui_div(rf, 1, xf, MPFR_RNDZ);

    NFLOAT_EXP(res) = -NFLOAT_EXP(x) + rf->_mpfr_exp;
    NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
nfloat_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    mpfr_t rf, xf, yf;
    slong prec = NFLOAT_CTX_PREC(ctx);

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_ZERO(x) && !NFLOAT_IS_SPECIAL(y))
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx); /* todo */
    }

    /* todo: make sure aliasing is correct */
    rf->_mpfr_d = NFLOAT_D(res);
    rf->_mpfr_prec = prec;
    rf->_mpfr_sign = 1;
    rf->_mpfr_exp = 0;

    xf->_mpfr_d = NFLOAT_D(x);
    xf->_mpfr_prec = prec;
    xf->_mpfr_sign = 1;
    xf->_mpfr_exp = 0;

    yf->_mpfr_d = NFLOAT_D(y);
    yf->_mpfr_prec = prec;
    yf->_mpfr_sign = 1;
    yf->_mpfr_exp = 0;

    mpfr_div(rf, xf, yf, MPFR_RNDZ);

    NFLOAT_EXP(res) = NFLOAT_EXP(x) - NFLOAT_EXP(y) + rf->_mpfr_exp;
    NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x) ^ NFLOAT_SGNBIT(y);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

/* check if mpfr_div_ui argument is ulong */
#if ULONG_MAX == UWORD_MAX

int
nfloat_div_ui(nfloat_ptr res, nfloat_srcptr x, ulong y, gr_ctx_t ctx)
{
    mpfr_t rf, xf;
    slong prec = NFLOAT_CTX_PREC(ctx);

    if (NFLOAT_IS_SPECIAL(x) || y == 0)
    {
        if (NFLOAT_IS_ZERO(x) && y != 0)
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx); /* todo */
    }

    /* todo: make sure aliasing is correct */
    rf->_mpfr_d = NFLOAT_D(res);
    rf->_mpfr_prec = prec;
    rf->_mpfr_sign = 1;
    rf->_mpfr_exp = 0;

    xf->_mpfr_d = NFLOAT_D(x);
    xf->_mpfr_prec = prec;
    xf->_mpfr_sign = 1;
    xf->_mpfr_exp = 0;

    mpfr_div_ui(rf, xf, y, MPFR_RNDZ);

    NFLOAT_EXP(res) = NFLOAT_EXP(x) + rf->_mpfr_exp;
    NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
nfloat_div_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx)
{
    mpfr_t rf, xf;
    slong prec = NFLOAT_CTX_PREC(ctx);

    if (NFLOAT_IS_SPECIAL(x) || y == 0)
    {
        if (NFLOAT_IS_ZERO(x) && y != 0)
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx); /* todo */
    }

    /* todo: make sure aliasing is correct */
    rf->_mpfr_d = NFLOAT_D(res);
    rf->_mpfr_prec = prec;
    rf->_mpfr_sign = 1;
    rf->_mpfr_exp = 0;

    xf->_mpfr_d = NFLOAT_D(x);
    xf->_mpfr_prec = prec;
    xf->_mpfr_sign = 1;
    xf->_mpfr_exp = 0;

    mpfr_div_si(rf, xf, y, MPFR_RNDZ);

    NFLOAT_EXP(res) = NFLOAT_EXP(x) + rf->_mpfr_exp;
    NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x) ^ (y < 0);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

#else

int
nfloat_div_ui(nfloat_ptr res, nfloat_srcptr x, ulong y, gr_ctx_t ctx)
{
    arf_t t;
    int status;

    arf_init(t);
    nfloat_get_arf(t, x, ctx);
    arf_div_ui(t, t, y, NFLOAT_CTX_PREC(ctx), ARF_RND_DOWN);
    status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return status;
}

int
nfloat_div_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx)
{
    arf_t t;
    int status;

    arf_init(t);
    nfloat_get_arf(t, x, ctx);
    arf_div_si(t, t, y, NFLOAT_CTX_PREC(ctx), ARF_RND_DOWN);
    status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return status;
}

#endif

int
nfloat_sqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    mpfr_t xf, zf;
    int odd_exp;
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);
    slong prec = NFLOAT_CTX_PREC(ctx);

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_nan(res, ctx);
        else
            return nfloat_set(res, x, ctx);
    }

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    odd_exp = NFLOAT_EXP(x) & 1;

    /* Powers of two */
    if (odd_exp &&
        NFLOAT_D(x)[nlimbs - 1] == (UWORD(1) << (FLINT_BITS - 1)) &&
        flint_mpn_zero_p(NFLOAT_D(x), nlimbs - 1))
    {
        nfloat_set(res, x, ctx);
        NFLOAT_EXP(res) = (NFLOAT_EXP(res) + 1) / 2;
        return GR_SUCCESS;
    }

    if (res == x)
    {
        zf->_mpfr_d = NFLOAT_D(x);
        zf->_mpfr_prec = prec;
        zf->_mpfr_sign = 1;
        zf->_mpfr_exp = odd_exp;

        mpfr_sqrt(zf, zf, MPFR_RNDZ);
    }
    else
    {
        xf->_mpfr_d = NFLOAT_D(x);
        xf->_mpfr_prec = prec;
        xf->_mpfr_sign = 1;
        xf->_mpfr_exp = odd_exp;

        zf->_mpfr_d = NFLOAT_D(res);
        zf->_mpfr_prec = prec;
        zf->_mpfr_sign = 1;
        zf->_mpfr_exp = 0;

        mpfr_sqrt(zf, xf, MPFR_RNDZ);
    }

    /* floor division */
    NFLOAT_EXP(res) = (NFLOAT_EXP(x) - (NFLOAT_EXP(x) < 0)) / 2 + zf->_mpfr_exp;

    NFLOAT_SGNBIT(res) = 0;

    return GR_SUCCESS;
}

int
nfloat_rsqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    mpfr_t xf, zf;
    int odd_exp;
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);
    slong prec = NFLOAT_CTX_PREC(ctx);

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_pos_inf(res, ctx);
        else if (NFLOAT_IS_POS_INF(x))
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx);
    }

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    odd_exp = NFLOAT_EXP(x) & 1;

    /* Powers of two */
    if (odd_exp &&
        NFLOAT_D(x)[nlimbs - 1] == (UWORD(1) << (FLINT_BITS - 1)) &&
        flint_mpn_zero_p(NFLOAT_D(x), nlimbs - 1))
    {
        nfloat_set(res, x, ctx);
        NFLOAT_EXP(res) = ((-NFLOAT_EXP(res) + 1)) / 2 + 1;
        return GR_SUCCESS;
    }

    if (res == x)
    {
        zf->_mpfr_d = NFLOAT_D(x);
        zf->_mpfr_prec = prec;
        zf->_mpfr_sign = 1;
        zf->_mpfr_exp = odd_exp;

        mpfr_rec_sqrt(zf, zf, MPFR_RNDZ);
    }
    else
    {
        xf->_mpfr_d = NFLOAT_D(x);
        xf->_mpfr_prec = prec;
        xf->_mpfr_sign = 1;
        xf->_mpfr_exp = odd_exp;

        zf->_mpfr_d = NFLOAT_D(res);
        zf->_mpfr_prec = prec;
        zf->_mpfr_sign = 1;
        zf->_mpfr_exp = 0;

        mpfr_rec_sqrt(zf, xf, MPFR_RNDZ);
    }

    /* floor division */
    NFLOAT_EXP(res) = -((NFLOAT_EXP(x) - (NFLOAT_EXP(x) < 0)) / 2) + zf->_mpfr_exp;

    NFLOAT_SGNBIT(res) = 0;

    return GR_SUCCESS;
}

static int
_nfloat_func1_via_arf(gr_method_unary_op func, nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    gr_ctx_t arf_ctx;
    arf_t t;
    int status;

    gr_ctx_init_real_float_arf(arf_ctx, NFLOAT_CTX_PREC(ctx));
    arf_init(t);
    nfloat_get_arf(t, x, ctx);
    status = func(t, t, arf_ctx);
    if (status == GR_SUCCESS)
        status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);
    gr_ctx_clear(arf_ctx);

    return status;
}

static int
_nfloat_func2_via_arf(gr_method_binary_op func, nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    gr_ctx_t arf_ctx;
    arf_t t, u;
    int status;

    gr_ctx_init_real_float_arf(arf_ctx, NFLOAT_CTX_PREC(ctx));
    arf_init(t);
    arf_init(u);
    nfloat_get_arf(t, x, ctx);
    nfloat_get_arf(u, y, ctx);
    status = func(t, t, u, arf_ctx);
    if (status == GR_SUCCESS)
        status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);
    arf_clear(u);
    gr_ctx_clear(arf_ctx);

    return status;
}

/* todo: fast code */
int nfloat_floor(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_floor, res, x, ctx); }
int nfloat_ceil(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_ceil, res, x, ctx); }
int nfloat_trunc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_trunc, res, x, ctx); }
int nfloat_nint(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_nint, res, x, ctx); }

int nfloat_pi(nfloat_ptr res, gr_ctx_t ctx)
{
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    FLINT_ASSERT(nlimbs <= ARB_PI4_TAB_LIMBS);

    NFLOAT_EXP(res) = 2;
    NFLOAT_SGNBIT(res) = 0;
    flint_mpn_copyi(NFLOAT_D(res), arb_pi4_tab + ARB_PI4_TAB_LIMBS - nlimbs, nlimbs);

    return GR_SUCCESS;
}

int nfloat_pow(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx) { return _nfloat_func2_via_arf((gr_method_binary_op) gr_pow, res, x, y, ctx); }
int nfloat_exp(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_exp, res, x, ctx); }
int nfloat_expm1(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_expm1, res, x, ctx); }
int nfloat_log(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_log, res, x, ctx); }
int nfloat_log1p(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_log1p, res, x, ctx); }
int nfloat_sin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_sin, res, x, ctx); }
int nfloat_cos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_cos, res, x, ctx); }
int nfloat_tan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_tan, res, x, ctx); }
int nfloat_sinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_sinh, res, x, ctx); }
int nfloat_cosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_cosh, res, x, ctx); }
int nfloat_tanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_tanh, res, x, ctx); }
int nfloat_atan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_atan, res, x, ctx); }
int nfloat_gamma(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_gamma, res, x, ctx); }
int nfloat_zeta(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_zeta, res, x, ctx); }
