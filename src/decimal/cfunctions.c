/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Elementary and special functions of complex decimal floats and balls.

    Real arguments are delegated to the real functions (which handle exact
    values, tiny and large arguments, and give an exact zero imaginary
    part). Purely imaginary arguments of the elementary functions are
    reduced to real functions through identities such as sin(iy) = i sinh(y).
    Tiny nonreal arguments are handled with Taylor expansions and a
    computed tail sign. Everything else goes through acb: directly for
    balls, and inside Ziv's loop, componentwise, for floats.
*/

#include <math.h>
#include "decimal.h"
#include "mag.h"
#include "fmpq.h"
#include "arf.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_special.h"

/* Not performance-critical: optimize for size. */
PUSH_OPTIONS
OPTIMIZE_OSIZE

#define RND(ctx) DECIMAL_CTX_RND(ctx)
#define RND_IM(ctx) DECIMAL_CTX_RND_IM(ctx)
#define PREC(ctx) DECIMAL_CTX_PREC(ctx)
#define EXACT_RND (DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS)

/* ------------------------------------------------------------------------- */
/*    Calling acb methods                                                    */
/* ------------------------------------------------------------------------- */

typedef struct
{
    int type;       /* 0: none, 1: ulong, 2: fmpz, 3: fmpq */
    ulong u;
    const fmpz * z;
    const fmpq * q;
}
extra_arg;

DECIMAL_DRIVER int
_gr_dispatch(gr_ptr * res, slong nres, gr_srcptr * args, slong nargs, int has_flag, int flag, const extra_arg * extra, gr_funcptr fn, gr_ctx_t actx)
{
    if (extra != NULL && extra->type != 0)
    {
        switch (nargs * 10 + extra->type)
        {
            case 1: return ((gr_method_unary_op_ui) fn)(res[0], extra->u, actx);
            case 2: return ((gr_method_unary_op_fmpz) fn)(res[0], extra->z, actx);
            case 3: return ((gr_method_unary_op_fmpq) fn)(res[0], extra->q, actx);
            case 11: return ((gr_method_binary_op_ui) fn)(res[0], args[0], extra->u, actx);
            case 12: return ((gr_method_binary_op_fmpz) fn)(res[0], args[0], extra->z, actx);
            case 13: return ((gr_method_binary_op_fmpq) fn)(res[0], args[0], extra->q, actx);
            default:
                flint_throw(FLINT_ERROR, "_gr_dispatch: unsupported signature\n");
        }
    }

    switch (nres * 100 + nargs * 10 + has_flag)
    {
        case 100: return ((gr_method_constant_op) fn)(res[0], actx);
        case 110: return ((gr_method_unary_op) fn)(res[0], args[0], actx);
        case 111: return ((gr_method_unary_op_with_flag) fn)(res[0], args[0], flag, actx);
        case 120: return ((gr_method_binary_op) fn)(res[0], args[0], args[1], actx);
        case 121: return ((gr_method_binary_op_with_flag) fn)(res[0], args[0], args[1], flag, actx);
        case 130: return ((gr_method_ternary_op) fn)(res[0], args[0], args[1], args[2], actx);
        case 131: return ((gr_method_ternary_op_with_flag) fn)(res[0], args[0], args[1], args[2], flag, actx);
        case 140: return ((gr_method_quaternary_op) fn)(res[0], args[0], args[1], args[2], args[3], actx);
        case 141: return ((gr_method_quaternary_op_with_flag) fn)(res[0], args[0], args[1], args[2], args[3], flag, actx);
        case 210: return ((gr_method_binary_unary_op) fn)(res[0], res[1], args[0], actx);
        case 211: return ((gr_method_binary_unary_op_with_flag) fn)(res[0], res[1], args[0], flag, actx);
        case 220: return ((gr_method_binary_binary_op) fn)(res[0], res[1], args[0], args[1], actx);
        case 410: return ((gr_method_quaternary_unary_op) fn)(res[0], res[1], res[2], res[3], args[0], actx);
        default:
            flint_throw(FLINT_ERROR, "_gr_dispatch: unsupported signature\n");
    }
}

#define MAX_ARGS 4
#define MAX_RES 4

DECIMAL_DRIVER int
_deccball_acb_gr_extra(deccball_ptr * res, slong nres, deccball_srcptr * args, slong nargs, int has_flag, int flag, const extra_arg * extra, int method, gr_ctx_t ctx)
{
    gr_ctx_t actx;
    acb_struct a[MAX_ARGS], r[MAX_RES];
    gr_ptr rp[MAX_RES];
    gr_srcptr asp[MAX_ARGS];
    decfloat_srcptr re[MAX_ARGS], im[MAX_ARGS];
    slong i, wp;
    int status = GR_SUCCESS;

    for (i = 0; i < nargs; i++)
    {
        re[i] = &args[i]->re.mid;
        im[i] = &args[i]->im.mid;
    }
    if (_decimal_arb_args_infeasible(method, re, im, nargs, PREC(ctx), ctx))
        return GR_UNABLE;

    wp = _decimal_digits_to_bits(PREC(ctx));
    gr_ctx_init_complex_acb(actx, wp);

    for (i = 0; i < nargs; i++)
    {
        acb_init(a + i);
        asp[i] = a + i;
        status |= deccball_get_acb(a + i, args[i], wp, ctx);
    }
    for (i = 0; i < nres; i++)
    {
        acb_init(r + i);
        rp[i] = r + i;
    }

    if (status == GR_SUCCESS)
        status = _gr_dispatch(rp, nres, asp, nargs, has_flag, flag, extra, actx->methods[method], actx);

    if (status == GR_SUCCESS)
        for (i = 0; i < nres; i++)
            status |= deccball_set_acb(res[i], r + i, ctx);

    for (i = 0; i < nargs; i++)
        acb_clear(a + i);
    for (i = 0; i < nres; i++)
        acb_clear(r + i);
    gr_ctx_clear(actx);
    return status;
}

/* rounds an acb result to a deccfloat; returns 1 if decided, 0 to retry
   at higher precision, -1 on error */
/* mask: bit 0 = real part, bit 1 = imaginary part; unmasked parts of
   res are left untouched */
static int
_deccfloat_round_acb(deccfloat_t res, const acb_t r, slong prec, int mask, gr_ctx_t ctx)
{
    decfloat_t tre, tim;
    int rr = 1;

    decfloat_init(tre, ctx);
    decfloat_init(tim, ctx);

    if (mask & 1)
        rr = _decfloat_round_arb(tre, acb_realref(r), prec, RND(ctx), ctx);
    if (rr == 1 && (mask & 2))
        rr = _decfloat_round_arb(tim, acb_imagref(r), prec, RND_IM(ctx), ctx);
    if (rr == 1)
    {
        if (mask & 1)
        {
            decfloat_swap(&res->re, tre, ctx);
            rr = (_decfloat_finalize(&res->re, ctx) == GR_SUCCESS) ? 1 : -1;
        }
        if (rr == 1 && (mask & 2))
        {
            decfloat_swap(&res->im, tim, ctx);
            rr = (_decfloat_finalize(&res->im, ctx) == GR_SUCCESS) ? 1 : -1;
        }
    }

    decfloat_clear(tre, ctx);
    decfloat_clear(tim, ctx);
    return rr;
}

/* acb method on floats, correctly rounded with Ziv's strategy */
/* mask (for a single result): which parts still need to be computed;
   the other parts of res[0] are kept */
DECIMAL_DRIVER int
_deccfloat_acb_gr_masked(deccfloat_ptr * res, slong nres, deccfloat_srcptr * args, slong nargs, int has_flag, int flag, const extra_arg * extra, int method, int mask, gr_ctx_t ctx)
{
    slong prec = PREC(ctx);
    gr_ctx_t actx;
    acb_struct a[MAX_ARGS], r[MAX_RES];
    gr_ptr rp[MAX_RES];
    gr_srcptr asp[MAX_ARGS];
    deccfloat_struct tmp[MAX_RES];
    slong i, wp, wpbits, wp_max;
    int status = GR_UNABLE, rr;

    for (i = 0; i < nargs; i++)
    {
        if (_deccfloat_is_nan(args[i]))
        {
            for (i = 0; i < nres; i++)
                status = deccfloat_nan(res[i], ctx);
            return status;
        }
    }

    if (prec == DECIMAL_PREC_EXACT)
        return GR_UNABLE;

    {
        decfloat_srcptr re[MAX_ARGS], im[MAX_ARGS];
        for (i = 0; i < nargs; i++)
        {
            re[i] = &args[i]->re;
            im[i] = &args[i]->im;
        }
        if (_decimal_arb_args_infeasible(method, re, im, nargs, prec, ctx))
            return GR_UNABLE;
    }

    for (i = 0; i < nargs; i++)
    {
        acb_init(a + i);
        asp[i] = a + i;
    }
    for (i = 0; i < nres; i++)
    {
        acb_init(r + i);
        rp[i] = r + i;
        deccfloat_init(tmp + i, ctx);
    }

    /* parts already computed by the caller */
    if (mask != 3)
    {
        if (!(mask & 1))
            GR_MUST_SUCCEED(decfloat_set_round(&tmp[0].re, &res[0]->re, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
        if (!(mask & 2))
            GR_MUST_SUCCEED(decfloat_set_round(&tmp[0].im, &res[0]->im, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
    }

    wp_max = 0;
    if (_decimal_method_is_elementary(method))
        for (i = 0; i < nargs; i++)
            wp_max = FLINT_MAX(wp_max, FLINT_MAX(decfloat_digits(&args[i]->re, ctx), decfloat_digits(&args[i]->im, ctx)));
    wp_max = _decimal_ziv_wp_max(prec, wp_max);

    for (wp = prec + 10; wp <= wp_max; wp *= 2)
    {
        wpbits = _decimal_digits_to_bits(wp);
        gr_ctx_init_complex_acb(actx, wpbits);

        status = GR_SUCCESS;
        for (i = 0; i < nargs; i++)
            status |= deccfloat_get_acb(a + i, args[i], wpbits, ctx);

        if (status == GR_SUCCESS)
            status = _gr_dispatch(rp, nres, asp, nargs, has_flag, flag, extra, actx->methods[method], actx);

        gr_ctx_clear(actx);

        if (status != GR_SUCCESS)
        {
            if (status & GR_DOMAIN)
                break;
            continue;
        }

        rr = 1;
        for (i = 0; i < nres && rr == 1; i++)
            rr = _deccfloat_round_acb(tmp + i, r + i, prec, (i == 0) ? mask : 3, ctx);

        if (rr == 1)
        {
            for (i = 0; i < nres; i++)
                deccfloat_swap(res[i], tmp + i, ctx);
            status = GR_SUCCESS;
            break;
        }

        status = GR_UNABLE;
        if (rr == -1)
            break;
    }

    for (i = 0; i < nargs; i++)
        acb_clear(a + i);
    for (i = 0; i < nres; i++)
    {
        acb_clear(r + i);
        deccfloat_clear(tmp + i, ctx);
    }
    return status;
}

int
_deccball_acb_gr(deccball_ptr * res, slong nres, deccball_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx)
{
    return _deccball_acb_gr_extra(res, nres, args, nargs, has_flag, flag, NULL, method, ctx);
}

DECIMAL_DRIVER int
_deccfloat_acb_gr_extra(deccfloat_ptr * res, slong nres, deccfloat_srcptr * args, slong nargs, int has_flag, int flag, const extra_arg * extra, int method, gr_ctx_t ctx)
{
    return _deccfloat_acb_gr_masked(res, nres, args, nargs, has_flag, flag, extra, method, 3, ctx);
}

int
_deccfloat_acb_gr(deccfloat_ptr * res, slong nres, deccfloat_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx)
{
    return _deccfloat_acb_gr_masked(res, nres, args, nargs, has_flag, flag, NULL, method, 3, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Real and imaginary arguments                                          */
/* ------------------------------------------------------------------------- */

/* calls the real function given by method with a given rounding mode for its result */
static int
_real_method_rnd(decfloat_t res, int method, const decfloat_t x, int rnd, gr_ctx_t ctx)
{
    int saved = RND(ctx), status;
    RND(ctx) = rnd;
    status = _decfloat_unary_method(res, x, method, ctx);
    RND(ctx) = saved;
    return status;
}

/* the same, negating the result: rounding is mirrored */
static int
_real_method_rnd_neg(decfloat_t res, int method, const decfloat_t x, int rnd, gr_ctx_t ctx)
{
    int status = _real_method_rnd(res, method, x, DECIMAL_RND_NEGATE(rnd), ctx);
    if (status == GR_SUCCESS)
    {
        if (DECFLOAT_IS_POS_INF(res)) _decfloat_neg_inf(res);
        else if (DECFLOAT_IS_NEG_INF(res)) _decfloat_pos_inf(res);
        else res->m.size = -res->m.size;
    }
    return status;
}

/* what to do with real arguments */
enum
{
    RK_NONE,        /* no real reduction */
    RK_FN,          /* the real function */
    RK_LOG,         /* log(-t) = log(t) + pi i */
    RK_ACOS,        /* acos(t) = i acosh(t) (t > 1), pi - i acosh(-t) (t < -1) */
    RK_ACOSH,       /* acosh(t) = i acos(t) (-1 <= t < 1), acosh(-t) + pi i (t < -1) */
    RK_EXP_PI_I     /* exp(pi i t) = cos(pi t) + i sin(pi t) */
};

/* what to do with imaginary arguments: f(iy) = g(y) or i h(y) or -i h(y) */
enum
{
    IM_NONE,
    IM_EXP,     /* cos y + i sin y */
    IM_SIN,     /* i sinh y */
    IM_COS,     /* cosh y */
    IM_TAN,     /* i tanh y */
    IM_COT,     /* -i coth y */
    IM_SEC,     /* sech y */
    IM_CSC,     /* -i csch y */
    IM_SINH,    /* i sin y */
    IM_COSH,    /* cos y */
    IM_TANH,    /* i tan y */
    IM_COTH,    /* -i cot y */
    IM_SECH,    /* sec y */
    IM_CSCH,    /* -i csc y */
    IM_ASIN,    /* i asinh y */
    IM_ATAN,    /* i atanh y */
    IM_ASINH,   /* i asin y */
    IM_ATANH,   /* i atan y */
    IM_ERF,     /* i erfi y */
    IM_ERFI,    /* i erf y */
    IM_SI,      /* i Shi y */
    IM_SHI      /* i Si y */
};

/* returns GR_SUCCESS, GR_DOMAIN (not applicable), or GR_UNABLE */
DECIMAL_DRIVER int
_deccfloat_imag_arg(deccfloat_t res, const decfloat_t y, int kind, gr_ctx_t ctx)
{
    int fre = -1, fim = -1;
    int neg = 0, status;
    decfloat_t t, u;

    switch (kind)
    {
        case IM_EXP: fre = GR_METHOD_COS; fim = GR_METHOD_SIN; break;
        case IM_SIN: fim = GR_METHOD_SINH; break;
        case IM_COS: fre = GR_METHOD_COSH; break;
        case IM_TAN: fim = GR_METHOD_TANH; break;
        case IM_COT: fim = GR_METHOD_COTH; neg = 1; break;
        case IM_SEC: fre = GR_METHOD_SECH; break;
        case IM_CSC: fim = GR_METHOD_CSCH; neg = 1; break;
        case IM_SINH: fim = GR_METHOD_SIN; break;
        case IM_COSH: fre = GR_METHOD_COS; break;
        case IM_TANH: fim = GR_METHOD_TAN; break;
        case IM_COTH: fim = GR_METHOD_COT; neg = 1; break;
        case IM_SECH: fre = GR_METHOD_SEC; break;
        case IM_CSCH: fim = GR_METHOD_CSC; neg = 1; break;
        case IM_ASIN: fim = GR_METHOD_ASINH; break;
        case IM_ATAN: fim = GR_METHOD_ATANH; break;
        case IM_ASINH: fim = GR_METHOD_ASIN; break;
        case IM_ATANH: fim = GR_METHOD_ATAN; break;
        case IM_ERF: fim = GR_METHOD_ERFI; break;
        case IM_ERFI: fim = GR_METHOD_ERF; break;
        case IM_SI: fim = GR_METHOD_SINH_INTEGRAL; break;
        case IM_SHI: fim = GR_METHOD_SIN_INTEGRAL; break;
        default: return GR_DOMAIN;
    }

    decfloat_init(t, ctx);
    decfloat_init(u, ctx);

    status = GR_SUCCESS;
    if (fre >= 0)
        status = _real_method_rnd(t, fre, y, RND(ctx), ctx);
    if (status == GR_SUCCESS && fim >= 0)
    {
        if (neg)
            status = _real_method_rnd_neg(u, fim, y, RND_IM(ctx), ctx);
        else
            status = _real_method_rnd(u, fim, y, RND_IM(ctx), ctx);
    }

    if (status == GR_SUCCESS)
    {
        decfloat_swap(&res->re, t, ctx);
        decfloat_swap(&res->im, u, ctx);
    }

    decfloat_clear(t, ctx);
    decfloat_clear(u, ctx);
    return status;
}

/* pi rounded in a given mode */
static int
_pi_rnd(decfloat_t res, int rnd, gr_ctx_t ctx)
{
    int saved = RND(ctx), status;
    RND(ctx) = rnd;
    status = decfloat_pi(res, ctx);
    RND(ctx) = saved;
    return status;
}

/* real argument: returns GR_SUCCESS, GR_DOMAIN (use acb), or GR_UNABLE */
DECIMAL_DRIVER int
_deccfloat_real_arg(deccfloat_t res, const decfloat_t x, int method, int rkind, gr_ctx_t ctx)
{
    decfloat_t t, u;
    int status = GR_DOMAIN;

    if (rkind == RK_NONE)
        return GR_DOMAIN;

    decfloat_init(t, ctx);
    decfloat_init(u, ctx);

    if (rkind == RK_EXP_PI_I)
    {
        status = _real_method_rnd(t, GR_METHOD_COS_PI, x, RND(ctx), ctx);
        if (status == GR_SUCCESS)
            status = _real_method_rnd(u, GR_METHOD_SIN_PI, x, RND_IM(ctx), ctx);
        if (status == GR_SUCCESS)
        {
            decfloat_swap(&res->re, t, ctx);
            decfloat_swap(&res->im, u, ctx);
        }
        goto cleanup;
    }

    /* complex values on the real line */
    if (DECFLOAT_IS_FINITE(x))
    {
        if (rkind == RK_LOG && _decfloat_sgn(x, ctx) < 0)
        {
            /* log(-t) = log(t) + pi i */
            GR_MUST_SUCCEED(decfloat_neg_round(u, x, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
            status = _real_method_rnd(t, GR_METHOD_LOG, u, RND(ctx), ctx);
            if (status == GR_SUCCESS)
                status = _pi_rnd(u, RND_IM(ctx), ctx);
            if (status == GR_SUCCESS)
            {
                decfloat_swap(&res->re, t, ctx);
                decfloat_swap(&res->im, u, ctx);
            }
            goto cleanup;
        }

        if (rkind == RK_ACOS && _decfloat_cmp_si(x, 1, ctx) > 0)
        {
            /* i acosh(t) */
            status = _real_method_rnd(u, GR_METHOD_ACOSH, x, RND_IM(ctx), ctx);
            if (status == GR_SUCCESS)
            {
                decfloat_zero(&res->re, ctx);
                decfloat_swap(&res->im, u, ctx);
            }
            goto cleanup;
        }

        if (rkind == RK_ACOS && _decfloat_cmp_si(x, -1, ctx) < 0)
        {
            /* pi - i acosh(-t) */
            GR_MUST_SUCCEED(decfloat_neg_round(u, x, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
            status = _real_method_rnd_neg(u, GR_METHOD_ACOSH, u, RND_IM(ctx), ctx);
            if (status == GR_SUCCESS)
                status = _pi_rnd(t, RND(ctx), ctx);
            if (status == GR_SUCCESS)
            {
                decfloat_swap(&res->re, t, ctx);
                decfloat_swap(&res->im, u, ctx);
            }
            goto cleanup;
        }

        if (rkind == RK_ACOSH && _decfloat_cmp_si(x, 1, ctx) < 0)
        {
            if (_decfloat_cmp_si(x, -1, ctx) >= 0)
            {
                /* i acos(t) */
                status = _real_method_rnd(u, GR_METHOD_ACOS, x, RND_IM(ctx), ctx);
                if (status == GR_SUCCESS)
                {
                    decfloat_zero(&res->re, ctx);
                    decfloat_swap(&res->im, u, ctx);
                }
            }
            else
            {
                /* acosh(-t) + pi i */
                GR_MUST_SUCCEED(decfloat_neg_round(u, x, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
                status = _real_method_rnd(t, GR_METHOD_ACOSH, u, RND(ctx), ctx);
                if (status == GR_SUCCESS)
                    status = _pi_rnd(u, RND_IM(ctx), ctx);
                if (status == GR_SUCCESS)
                {
                    decfloat_swap(&res->re, t, ctx);
                    decfloat_swap(&res->im, u, ctx);
                }
            }
            goto cleanup;
        }
    }

    status = _decfloat_unary_method(t, x, method, ctx);

    if (status == GR_SUCCESS)
    {
        decfloat_swap(&res->re, t, ctx);
        decfloat_zero(&res->im, ctx);
    }
    else if (status == GR_DOMAIN && DECFLOAT_IS_ZERO(x) && rkind == RK_LOG)
    {
        /* a pole: no point in trying acb */
        status = GR_DOMAIN | GR_UNABLE;
    }

cleanup:
    decfloat_clear(t, ctx);
    decfloat_clear(u, ctx);
    return status;
}

/* ------------------------------------------------------------------------- */
/*    Tiny arguments                                                         */
/* ------------------------------------------------------------------------- */

/*
    f(z) = S + t where S is an exact leading part (z, 1, +/- 1/z, z - 1 or
    1/(z - 1)) and t is a tail bounded by 10^(tmul E + tadd), with
    |z| < sqrt(2) 10^E (E = the larger component exponent plus one), so
    that |z|^k < 10^(kE + 1) for k <= 6.

    The sign of each component of t is found from the leading terms of the
    Taylor expansion t = sum_j c_j z^(k_j) + ...: the first term with a
    nonzero component decides, provided that it dominates the bound
    (rnum_j / rden_j) |z|^(rk_j) on everything after it. Coefficients are
    given by their sign and a lower bound |c_j| >= clo_num / clo_den (exact
    when the coefficient is rational, which allows a term to serve as the
    leading part of a component whose S is zero). The powers of z are
    computed in ball arithmetic.
*/

enum { S_X, S_ONE, S_INVX, S_NEGINVX, S_XM1, S_INVXM1 };

typedef struct
{
    int csign;
    slong clo_num;
    slong clo_den;
    int exact;      /* coefficient is exactly clo_num / clo_den */
    int k;          /* power of z */
    slong rnum, rden;   /* the sum of all later terms is bounded by */
    int rk;             /* (rnum / rden) |z|^rk */
}
cterm;

typedef struct
{
    int kind;
    int tmul, tadd;
    int nterms;
    cterm terms[3];
}
ctail_info;

/* rounds (num/den) w + t to prec digits, where w is exact and t has a
   known sign with |t| < 10^tail_exp; returns 1, 0 or -1 */
static int
_round_scaled_with_tail(decfloat_t res, const decfloat_t w, slong num, slong den, int tail_sign, slong tail_exp, slong prec, int rnd, gr_ctx_t ctx)
{
    decfloat_t T, D, S;
    decimal_rounding_info info;
    slong P, ES;
    int r, status;

    decfloat_init(T, ctx);
    decfloat_init(D, ctx);
    decfloat_init(S, ctx);

    GR_MUST_SUCCEED(decfloat_set_round_si(T, num, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
    GR_MUST_SUCCEED(decfloat_set_round_si(D, den, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
    status = _decfloat_mul(T, T, w, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);

    /* S = T/D rounded toward -tail_sign at a high precision, so that the
       division error has the sign of the tail */
    P = prec + 4 * DECIMAL_CTX_E(ctx) + 10;
    if (status == GR_SUCCESS)
        status = _decfloat_div(S, T, D, P, ((tail_sign > 0) ? DECIMAL_RND_FLOOR : DECIMAL_RND_CEIL) | DECIMAL_RND_NOLIMITS, &info, NULL, ctx);

    if (status != GR_SUCCESS || DECFLOAT_IS_SPECIAL(S))
    {
        r = 0;
    }
    else
    {
        if (info.inexact)
        {
            ES = _decfloat_sci_exp_clamped(S, ctx);
            tail_exp = FLINT_MAX(tail_exp, ES - P + 2) + 1;
        }
        r = _decfloat_round_with_tail(res, S, tail_sign, tail_exp, prec, rnd, ctx);
    }

    decfloat_clear(T, ctx);
    decfloat_clear(D, ctx);
    decfloat_clear(S, ctx);
    return r;
}

/* whether (num/den) lb > B */
static int
_lb_dominates(const decmag_t lb, slong num, slong den, const decmag_t B, gr_ctx_t bctx)
{
    decmag_t L, R, t;
    int result;
    _decmag_init(L, bctx);
    _decmag_init(R, bctx);
    _decmag_init(t, bctx);
    _decmag_set_ui_lower(t, num, bctx);
    _decmag_mul_lower(L, lb, t, bctx);
    _decmag_mul_ui(R, B, den, bctx);
    result = _decmag_cmp(L, R, bctx) > 0;
    _decmag_clear(L, bctx);
    _decmag_clear(R, bctx);
    _decmag_clear(t, bctx);
    return result;
}

/* an exponent e with B < 10^e */
static slong
_decmag_exp_bound(const decmag_t B, gr_ctx_t bctx)
{
    slong e;
    if (DECMAG_IS_ZERO(B))
        return -DECFLOAT_EXP_CLAMP;
    if (DECMAG_IS_INF(B))
        return DECFLOAT_EXP_CLAMP;
    /* B = m 10^exp < 10^(exp + digits) */
    if (fmpz_fits_si(&B->exp))
    {
        e = fmpz_get_si(&B->exp);
        e = FLINT_MAX(FLINT_MIN(e, DECFLOAT_EXP_CLAMP / 2), -DECFLOAT_EXP_CLAMP / 2);
        return e + _decmag_digits(B);
    }
    return (fmpz_sgn(&B->exp) > 0) ? DECFLOAT_EXP_CLAMP : -DECFLOAT_EXP_CLAMP;
}

/* sign of a component of a power computed in ball arithmetic, with a
   lower bound for its absolute value: returns the sign (0 for an exact
   zero, or 2 if the sign is undetermined) */
static int
_ball_sign_lbound(decmag_t lb, const decball_t x, gr_ctx_t bctx)
{
    if (DECFLOAT_IS_ZERO(&x->mid) && DECMAG_IS_ZERO(&x->rad))
        return 0;
    if (!_decball_is_finite(x, bctx) || _decball_contains_zero(x, bctx))
        return 2;
    decball_get_abs_lbound(lb, x, bctx);
    return _decfloat_sgn(&x->mid, bctx);
}

/* returns a mask of the handled parts (bit 0: real, bit 1: imaginary),
   which are written to res, or -1 on error */
DECIMAL_DRIVER int
_deccfloat_tiny_arg(deccfloat_t res, const deccfloat_t z, const ctail_info * ti, slong prec, int rnd, int rnd_im, gr_ctx_t ctx)
{
    gr_ctx_t bctx;
    deccfloat_t h, S;
    deccball_t B, P[8];
    decfloat_t tre, tim;
    decmag_t lb, ub, ubt;
    const deccfloat_struct * base;
    slong E, i, npow = 0;
    int r = 0, status, comp;

    if (!_deccfloat_is_finite(z) || prec == DECIMAL_PREC_EXACT)
        return 0;

    if (DECFLOAT_IS_ZERO(&z->re) && DECFLOAT_IS_ZERO(&z->im))
        return 0;

    /* ball context for the powers: exact when the results are
       representable, otherwise with a rigorous radius */
    _gr_ctx_init_decimal(bctx, DECIMAL_CTX_CBALL, DECIMAL_CTX_E(ctx), prec + 4 * DECIMAL_CTX_E(ctx) + 40, DECIMAL_RND_NEAR, 0);
    decimal_ctx_set_rad_prec(bctx, DECMAG_MAX_PREC);

    deccfloat_init(h, ctx);
    deccfloat_init(S, ctx);
    deccball_init(B, bctx);
    decfloat_init(tre, ctx);
    decfloat_init(tim, ctx);
    _decmag_init(lb, bctx);
    _decmag_init(ub, bctx);
    _decmag_init(ubt, bctx);

    if (ti->kind == S_XM1 || ti->kind == S_INVXM1)
    {
        deccfloat_one(S, ctx);
        if (_deccfloat_add_exact(h, z, S, 1, ctx) != GR_SUCCESS)
            goto cleanup;
        base = h;
    }
    else
        base = z;

    if (DECFLOAT_IS_ZERO(&base->re) && DECFLOAT_IS_ZERO(&base->im))
        goto cleanup;

    E = WORD_MIN;
    if (!DECFLOAT_IS_ZERO(&base->re))
        E = FLINT_MAX(E, _decfloat_sci_exp_clamped(&base->re, ctx));
    if (!DECFLOAT_IS_ZERO(&base->im))
        E = FLINT_MAX(E, _decfloat_sci_exp_clamped(&base->im, ctx));
    E = E + 1;

    if (E > -2)
        goto cleanup;

    switch (ti->kind)
    {
        case S_X:
        case S_XM1:
            status = decfloat_set_round(&S->re, &base->re, DECIMAL_PREC_EXACT, EXACT_RND, ctx);
            status |= decfloat_set_round(&S->im, &base->im, DECIMAL_PREC_EXACT, EXACT_RND, ctx);
            break;
        case S_ONE:
            status = deccfloat_one(S, ctx);
            break;
        case S_INVX:
        case S_NEGINVX:
        case S_INVXM1:
            status = _deccfloat_inv_exact(S, base, ctx);
            if (status == GR_SUCCESS && ti->kind == S_NEGINVX)
            {
                S->re.m.size = -S->re.m.size;
                S->im.m.size = -S->im.m.size;
            }
            break;
        default:
            status = GR_UNABLE;
    }

    if (status != GR_SUCCESS)
        goto cleanup;

    /* upper bound for |z| */
    _decmag_set_decfloat(ub, &base->re, bctx);
    _decmag_set_decfloat(ubt, &base->im, bctx);
    _decmag_hypot(ub, ub, ubt, bctx);

    /* powers base^k for the terms */
    {
        int kmax = 0;
        for (i = 0; i < ti->nterms; i++)
            kmax = FLINT_MAX(kmax, ti->terms[i].k);
        if (deccball_set_deccfloat(B, base, bctx) != GR_SUCCESS)
            goto cleanup;
        for (i = 0; i <= kmax; i++)
        {
            deccball_init(P[i], bctx);
            npow = i + 1;
            if (i == 0)
                deccball_one(P[i], bctx);
            else if (deccball_mul(P[i], P[i - 1], B, bctx) != GR_SUCCESS)
                goto cleanup;
        }
    }

    for (comp = 0; comp < 2; comp++)
    {
        const decfloat_struct * Sc = comp ? &S->im : &S->re;
        decfloat_ptr out = comp ? tim : tre;
        int crnd = comp ? rnd_im : rnd;
        int sgn[3], found = 0, tsign = 0, rr;
        slong jfound = 0;
        decmag_struct lbc[3], ubc[3], after[3];
        slong j, n = ti->nterms;

        for (j = 0; j < n; j++)
        {
            const decball_struct * Pc = comp ? &P[ti->terms[j].k]->im : &P[ti->terms[j].k]->re;
            _decmag_init(lbc + j, bctx);
            _decmag_init(ubc + j, bctx);
            sgn[j] = _ball_sign_lbound(lbc + j, Pc, bctx);
            decball_get_abs_ubound(ubc + j, Pc, bctx);
        }

        /* after[j] bounds the sum of the terms after term j in this
           component: the remainder after the last term (from |z|) plus
           the listed terms j+1, ..., n-1 */
        for (j = 0; j < n; j++)
            _decmag_init(after + j, bctx);
        _decmag_pow_ui(after + n - 1, ub, ti->terms[n - 1].rk, bctx);
        _decmag_mul_ui(after + n - 1, after + n - 1, ti->terms[n - 1].rnum, bctx);
        _decmag_div_ui(after + n - 1, after + n - 1, ti->terms[n - 1].rden, bctx);
        for (j = n - 2; j >= 0; j--)
        {
            const cterm * T = ti->terms + j + 1;
            decmag_t c;
            _decmag_init(c, bctx);
            /* |c_{j+1}| <= num/den (exact) or 2 num/den (a lower bound was given) */
            _decmag_set_ui(c, T->exact ? T->clo_num : 2 * T->clo_num, bctx);
            _decmag_div_ui(c, c, T->clo_den, bctx);
            _decmag_mul(c, c, ubc + j + 1, bctx);
            _decmag_add(after + j, after + j + 1, c, bctx);
            _decmag_clear(c, bctx);
        }

        /* the first term with a nonzero component decides the sign, if
           it dominates everything after it */
        for (j = 0; j < n; j++)
        {
            const cterm * T = ti->terms + j;

            if (sgn[j] == 0)
                continue;
            if (sgn[j] == 2)
                break;

            if (DECFLOAT_IS_ZERO(Sc))
            {
                /* this term is the leading part (it must be exact); the
                   next term with a nonzero component gives the sign of the
                   remainder */
                const decball_struct * Pc = comp ? &P[T->k]->im : &P[T->k]->re;
                slong j2;
                int found2 = 0, tsign2 = 0;

                if (!T->exact || !DECMAG_IS_ZERO(&Pc->rad))
                    break;

                for (j2 = j + 1; j2 < n; j2++)
                {
                    if (sgn[j2] == 0)
                        continue;
                    if (sgn[j2] == 2)
                        break;
                    if (_lb_dominates(lbc + j2, ti->terms[j2].clo_num, ti->terms[j2].clo_den, after + j2, bctx))
                    {
                        tsign2 = ti->terms[j2].csign * sgn[j2];
                        found2 = 1;
                    }
                    break;
                }

                if (found2)
                {
                    rr = _round_scaled_with_tail(out, &Pc->mid, T->csign * T->clo_num, T->clo_den, tsign2, _decmag_exp_bound(after + j, bctx), prec, crnd, ctx);
                    if (rr == -1) { r = -1; goto cleanup2; }
                    found = (rr == 1) ? 2 : 0;
                }
                break;
            }

            if (_lb_dominates(lbc + j, T->clo_num, T->clo_den, after + j, bctx))
            {
                tsign = T->csign * sgn[j];
                found = 1;
                jfound = j;
            }
            break;
        }

        if (found == 1)
        {
            /* the tail is term jfound plus everything after it */
            const cterm * T = ti->terms + jfound;
            decmag_t tb;
            slong tail_exp;

            _decmag_init(tb, bctx);
            _decmag_set_ui(tb, T->exact ? T->clo_num : 2 * T->clo_num, bctx);
            _decmag_div_ui(tb, tb, T->clo_den, bctx);
            _decmag_mul(tb, tb, ubc + jfound, bctx);
            _decmag_add(tb, tb, after + jfound, bctx);
            tail_exp = _decmag_exp_bound(tb, bctx);
            _decmag_clear(tb, bctx);

            rr = _decfloat_round_with_tail(out, Sc, tsign, tail_exp, prec, crnd, ctx);
            if (rr == -1) { r = -1; goto cleanup2; }
            if (rr == 1) r |= (1 << comp);
        }
        else if (found == 2)
            r |= (1 << comp);

cleanup2:
        for (j = 0; j < n; j++)
        {
            _decmag_clear(lbc + j, bctx);
            _decmag_clear(ubc + j, bctx);
        }
        for (j = 0; j < n; j++)
            _decmag_clear(after + j, bctx);
        if (r == -1)
            goto cleanup;
    }

    if (r & 1)
        decfloat_swap(&res->re, tre, ctx);
    if (r & 2)
        decfloat_swap(&res->im, tim, ctx);

cleanup:
    for (i = 0; i < npow; i++)
        deccball_clear(P[i], bctx);
    deccfloat_clear(h, ctx);
    deccfloat_clear(S, ctx);
    deccball_clear(B, bctx);
    decfloat_clear(tre, ctx);
    decfloat_clear(tim, ctx);
    _decmag_clear(lb, bctx);
    _decmag_clear(ub, bctx);
    _decmag_clear(ubt, bctx);
    gr_ctx_clear(bctx);
    return r;
}

#define CT(name, kind, tmul, tadd, ...) static const ctail_info name = {kind, tmul, tadd, __VA_ARGS__};

/* tables: {kind, tmul, tadd, nterms, {csign, num, den, exact, k, rnum, rden, rk}...},
   where the remainder after a term is bounded by (rnum / rden) |z|^rk */
CT(ct_exp,     S_ONE,    1, 1, 2, {{ 1, 1, 1, 1, 1, 51, 100, 2}, { 1, 1, 2, 1, 2, 17, 100, 3}})
CT(ct_expm1,   S_X,      2, 1, 2, {{ 1, 1, 2, 1, 2, 17, 100, 3}, { 1, 1, 6, 1, 3, 5, 100, 4}})
CT(ct_log1p,   S_X,      2, 1, 2, {{-1, 1, 2, 1, 2, 34, 100, 3}, { 1, 1, 3, 1, 3, 26, 100, 4}})
CT(ct_log,     S_XM1,    2, 1, 2, {{-1, 1, 2, 1, 2, 34, 100, 3}, { 1, 1, 3, 1, 3, 26, 100, 4}})
CT(ct_sin,     S_X,      3, 0, 2, {{-1, 1, 6, 1, 3, 1, 100, 5}, { 1, 1, 120, 1, 5, 1, 1000, 7}})
CT(ct_cos,     S_ONE,    2, 1, 3, {{-1, 1, 2, 1, 2, 5, 100, 4}, { 1, 1, 24, 1, 4, 2, 1000, 6}, {-1, 1, 720, 1, 6, 3, 100000, 8}})
CT(ct_tan,     S_X,      3, 0, 2, {{ 1, 1, 3, 1, 3, 14, 100, 5}, { 1, 2, 15, 1, 5, 6, 100, 7}})
CT(ct_sinh,    S_X,      3, 0, 2, {{ 1, 1, 6, 1, 3, 1, 100, 5}, { 1, 1, 120, 1, 5, 1, 1000, 7}})
CT(ct_cosh,    S_ONE,    2, 1, 3, {{ 1, 1, 2, 1, 2, 5, 100, 4}, { 1, 1, 24, 1, 4, 2, 1000, 6}, { 1, 1, 720, 1, 6, 3, 100000, 8}})
CT(ct_tanh,    S_X,      3, 0, 2, {{-1, 1, 3, 1, 3, 14, 100, 5}, { 1, 2, 15, 1, 5, 6, 100, 7}})
CT(ct_asin,    S_X,      3, 0, 2, {{ 1, 1, 6, 1, 3, 8, 100, 5}, { 1, 3, 40, 1, 5, 5, 100, 7}})
CT(ct_atan,    S_X,      3, 0, 2, {{-1, 1, 3, 1, 3, 21, 100, 5}, { 1, 1, 5, 1, 5, 15, 100, 7}})
CT(ct_asinh,   S_X,      3, 0, 2, {{-1, 1, 6, 1, 3, 8, 100, 5}, { 1, 3, 40, 1, 5, 5, 100, 7}})
CT(ct_atanh,   S_X,      3, 0, 2, {{ 1, 1, 3, 1, 3, 21, 100, 5}, { 1, 1, 5, 1, 5, 15, 100, 7}})
CT(ct_sinc,    S_ONE,    2, 0, 3, {{-1, 1, 6, 1, 2, 1, 100, 4}, { 1, 1, 120, 1, 4, 2, 10000, 6}, {-1, 1, 5040, 1, 6, 3, 1000000, 8}})
CT(ct_sec,     S_ONE,    2, 1, 3, {{ 1, 1, 2, 1, 2, 22, 100, 4}, { 1, 5, 24, 1, 4, 9, 100, 6}, { 1, 61, 720, 1, 6, 4, 100, 8}})
CT(ct_sech,    S_ONE,    2, 1, 3, {{-1, 1, 2, 1, 2, 22, 100, 4}, { 1, 5, 24, 1, 4, 9, 100, 6}, {-1, 61, 720, 1, 6, 4, 100, 8}})
CT(ct_cot,     S_INVX,   1, 0, 2, {{-1, 1, 3, 1, 1, 23, 1000, 3}, {-1, 1, 45, 1, 3, 22, 10000, 5}})
CT(ct_csc,     S_INVX,   1, 0, 2, {{ 1, 1, 6, 1, 1, 2, 100, 3}, { 1, 7, 360, 1, 3, 21, 10000, 5}})
CT(ct_coth,    S_INVX,   1, 0, 2, {{ 1, 1, 3, 1, 1, 23, 1000, 3}, {-1, 1, 45, 1, 3, 22, 10000, 5}})
CT(ct_csch,    S_INVX,   1, 0, 2, {{-1, 1, 6, 1, 1, 2, 100, 3}, { 1, 7, 360, 1, 3, 21, 10000, 5}})
CT(ct_gamma,   S_INVX,   0, 0, 2, {{-1, 1, 2, 0, 0, 101, 100, 1}, { 1, 9, 10, 0, 1, 96, 100, 2}})
CT(ct_digamma, S_NEGINVX, 0, 0, 2, {{-1, 1, 2, 0, 0, 168, 100, 1}, { 1, 3, 2, 0, 1, 123, 100, 2}})
CT(ct_zeta,    S_INVXM1, 0, 0, 2, {{ 1, 1, 2, 0, 0, 8, 100, 1}, { 1, 7, 100, 0, 1, 6, 1000, 2}})
CT(ct_rgamma,  S_X,      2, 1, 2, {{ 1, 1, 2, 0, 2, 69, 100, 3}, {-1, 6, 10, 0, 3, 5, 100, 4}})
CT(ct_lambertw, S_X,     2, 1, 3, {{-1, 1, 1, 1, 2, 158, 100, 3}, { 1, 3, 2, 1, 3, 28, 10, 4}, {-1, 8, 3, 1, 4, 55, 10, 5}})
CT(ct_dilog,   S_X,      2, 0, 2, {{ 1, 1, 4, 1, 2, 12, 100, 3}, { 1, 1, 9, 1, 3, 7, 100, 4}})
CT(ct_cos_pi,  S_ONE,    2, 2, 2, {{-1, 4, 1, 0, 2, 41, 10, 4}, { 1, 4, 1, 0, 4, 14, 10, 6}})
CT(ct_sinc_pi, S_ONE,    2, 1, 2, {{-1, 3, 2, 0, 2, 82, 100, 4}, { 1, 4, 5, 0, 4, 2, 10, 6}})
CT(ct_si,      S_X,      3, 0, 2, {{-1, 1, 18, 1, 3, 2, 1000, 5}, { 1, 1, 600, 1, 5, 3, 100000, 7}})
CT(ct_shi,     S_X,      3, 0, 2, {{ 1, 1, 18, 1, 3, 2, 1000, 5}, { 1, 1, 600, 1, 5, 3, 100000, 7}})

/* ------------------------------------------------------------------------- */
/*    Unary functions                                                        */
/* ------------------------------------------------------------------------- */

DECIMAL_DRIVER int
_deccfloat_unary(deccfloat_t res, const deccfloat_t x, int method, int rkind, int ikind, const ctail_info * ti, gr_ctx_t ctx)
{
    int r;
    deccfloat_srcptr a[1];
    deccfloat_ptr rr[1];

    if (_deccfloat_is_nan(x))
        return deccfloat_nan(res, ctx);

    if (DECFLOAT_IS_ZERO(&x->im))
    {
        /* the real function may fail outside its real domain: try acb */
        r = _deccfloat_real_arg(res, &x->re, method, rkind, ctx);
        if (r == (GR_DOMAIN | GR_UNABLE))
            return GR_DOMAIN;
        if (r == GR_SUCCESS)
            return r;
    }
    else
    {
        if (DECFLOAT_IS_ZERO(&x->re) && ikind != IM_NONE)
        {
            /* the real function may fail outside its real domain: try acb */
            r = _deccfloat_imag_arg(res, &x->im, ikind, ctx);
            if (r == GR_SUCCESS)
                return r;
        }

        if (ti != NULL)
        {
            /* the parts not handled here are left to acb */
            r = _deccfloat_tiny_arg(res, x, ti, PREC(ctx), RND(ctx), RND_IM(ctx), ctx);
            if (r == 3) return GR_SUCCESS;
            if (r == -1) return GR_UNABLE;
            if (r != 0)
            {
                a[0] = x;
                rr[0] = res;
                return _deccfloat_acb_gr_masked(rr, 1, a, 1, 0, 0, NULL, method, 3 & ~r, ctx);
            }
        }
    }

    a[0] = x;
    rr[0] = res;
    return _deccfloat_acb_gr(rr, 1, a, 1, 0, 0, method, ctx);
}

DECIMAL_DRIVER int
_deccball_unary(deccball_t res, const deccball_t x, int method, int real, gr_ctx_t ctx)
{
    deccball_srcptr a[1];
    deccball_ptr rr[1];

    if (real && _deccball_is_real(x, ctx))
    {
        decball_t t;
        int status;
        decball_init(t, ctx);
        status = _decball_unary_fn(t, &x->re, method, ctx);
        if (status == GR_SUCCESS)
        {
            decball_swap(&res->re, t, ctx);
            decball_zero(&res->im, ctx);
        }
        decball_clear(t, ctx);
        if (status == GR_SUCCESS)
            return status;
    }

    a[0] = x;
    rr[0] = res;
    return _deccball_acb_gr(rr, 1, a, 1, 0, 0, method, ctx);
}

/* how each unary function treats real, imaginary, tiny and large arguments */
typedef struct
{
    int method;
    signed char rkind;          /* RK_* */
    signed char ikind;          /* IM_* */
    const ctail_info * tail;    /* tiny arguments (or tiny 1/z if inv) */
    signed char inv;            /* f(z) = g(1/z) for large z */
}
complex_spec;

static const complex_spec _complex_specs[] =
{
/*   method                             real         imaginary tail          inv */
    {GR_METHOD_EXP,                     RK_FN,       IM_EXP,   &ct_exp,      0},
    {GR_METHOD_EXPM1,                   RK_FN,       IM_NONE,  &ct_expm1,    0},
    {GR_METHOD_LOG,                     RK_LOG,      IM_NONE,  &ct_log,      0},
    {GR_METHOD_LOG1P,                   RK_FN,       IM_NONE,  &ct_log1p,    0},
    {GR_METHOD_SIN,                     RK_FN,       IM_SIN,   &ct_sin,      0},
    {GR_METHOD_COS,                     RK_FN,       IM_COS,   &ct_cos,      0},
    {GR_METHOD_TAN,                     RK_FN,       IM_TAN,   &ct_tan,      0},
    {GR_METHOD_COT,                     RK_FN,       IM_COT,   &ct_cot,      0},
    {GR_METHOD_SEC,                     RK_FN,       IM_SEC,   &ct_sec,      0},
    {GR_METHOD_CSC,                     RK_FN,       IM_CSC,   &ct_csc,      0},
    {GR_METHOD_SINC,                    RK_FN,       IM_NONE,  &ct_sinc,     0},
    {GR_METHOD_SIN_PI,                  RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_COS_PI,                  RK_FN,       IM_NONE,  &ct_cos_pi,   0},
    {GR_METHOD_TAN_PI,                  RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_COT_PI,                  RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_SEC_PI,                  RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_CSC_PI,                  RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_SINC_PI,                 RK_FN,       IM_NONE,  &ct_sinc_pi,  0},
    {GR_METHOD_SINH,                    RK_FN,       IM_SINH,  &ct_sinh,     0},
    {GR_METHOD_COSH,                    RK_FN,       IM_COSH,  &ct_cosh,     0},
    {GR_METHOD_TANH,                    RK_FN,       IM_TANH,  &ct_tanh,     0},
    {GR_METHOD_COTH,                    RK_FN,       IM_COTH,  &ct_coth,     0},
    {GR_METHOD_SECH,                    RK_FN,       IM_SECH,  &ct_sech,     0},
    {GR_METHOD_CSCH,                    RK_FN,       IM_CSCH,  &ct_csch,     0},
    {GR_METHOD_ASIN,                    RK_FN,       IM_ASIN,  &ct_asin,     0},
    {GR_METHOD_ACOS,                    RK_ACOS,     IM_NONE,  NULL,         0},
    {GR_METHOD_ATAN,                    RK_FN,       IM_ATAN,  &ct_atan,     0},
    {GR_METHOD_ASEC,                    RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ASINH,                   RK_FN,       IM_ASINH, &ct_asinh,    0},
    {GR_METHOD_ACOSH,                   RK_ACOSH,    IM_NONE,  NULL,         0},
    {GR_METHOD_ATANH,                   RK_FN,       IM_ATANH, &ct_atanh,    0},
    {GR_METHOD_ASECH,                   RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ASIN_PI,                 RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ACOS_PI,                 RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ATAN_PI,                 RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ACOT_PI,                 RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ASEC_PI,                 RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ACSC_PI,                 RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_LAMBERTW,                RK_FN,       IM_NONE,  &ct_lambertw, 0},
    {GR_METHOD_GAMMA,                   RK_FN,       IM_NONE,  &ct_gamma,    0},
    {GR_METHOD_RGAMMA,                  RK_FN,       IM_NONE,  &ct_rgamma,   0},
    {GR_METHOD_LGAMMA,                  RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_DIGAMMA,                 RK_FN,       IM_NONE,  &ct_digamma,  0},
    {GR_METHOD_ZETA,                    RK_FN,       IM_NONE,  &ct_zeta,     0},
    {GR_METHOD_BARNES_G,                RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_LOG_BARNES_G,            RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ERF,                     RK_FN,       IM_ERF,   NULL,         0},
    {GR_METHOD_ERFC,                    RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ERFI,                    RK_FN,       IM_ERFI,  NULL,         0},
    {GR_METHOD_ERFINV,                  RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ERFCINV,                 RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_EXP_INTEGRAL_EI,         RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_SIN_INTEGRAL,            RK_FN,       IM_SI,    &ct_si,       0},
    {GR_METHOD_COS_INTEGRAL,            RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_SINH_INTEGRAL,           RK_FN,       IM_SHI,   &ct_shi,      0},
    {GR_METHOD_COSH_INTEGRAL,           RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_DILOG,                   RK_FN,       IM_NONE,  &ct_dilog,    0},
    {GR_METHOD_AGM1,                    RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_AIRY_AI,                 RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_AIRY_BI,                 RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_AIRY_AI_PRIME,           RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_AIRY_BI_PRIME,           RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_LOG2,                    RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_LOG10,                   RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_EXP2,                    RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_EXP10,                   RK_FN,       IM_NONE,  NULL,         0},
    {GR_METHOD_ACOT,                    RK_FN,       IM_NONE,  &ct_atan,     1},
    {GR_METHOD_ACOTH,                   RK_FN,       IM_NONE,  &ct_atanh,    1},
    {GR_METHOD_ACSC,                    RK_FN,       IM_NONE,  &ct_asin,     1},
    {GR_METHOD_ACSCH,                   RK_FN,       IM_NONE,  &ct_asinh,    1},
    {GR_METHOD_EXP_PI_I,                RK_EXP_PI_I, IM_NONE,  NULL,         0},
    {GR_METHOD_LOG_PI_I,                RK_NONE,     IM_NONE,  NULL,         0},
    {GR_METHOD_DIRICHLET_ETA,           RK_NONE,     IM_NONE,  NULL,         0},
    {GR_METHOD_RIEMANN_XI,              RK_NONE,     IM_NONE,  NULL,         0},
    {GR_METHOD_MODULAR_J,               RK_NONE,     IM_NONE,  NULL,         0},
    {GR_METHOD_MODULAR_LAMBDA,          RK_NONE,     IM_NONE,  NULL,         0},
    {GR_METHOD_MODULAR_DELTA,           RK_NONE,     IM_NONE,  NULL,         0},
    {GR_METHOD_DEDEKIND_ETA,            RK_NONE,     IM_NONE,  NULL,         0},
    {GR_METHOD_ELLIPTIC_K,              RK_NONE,     IM_NONE,  NULL,         0},
    {GR_METHOD_ELLIPTIC_E,              RK_NONE,     IM_NONE,  NULL,         0},
    {-1, 0, 0, NULL, 0}
};

static const complex_spec *
_complex_spec(int method)
{
    static signed char index[GR_METHOD_TAB_SIZE];
    static volatile int initialized = 0;
    slong i;

    if (!initialized)
    {
        for (i = 0; i < GR_METHOD_TAB_SIZE; i++)
            index[i] = -1;
        for (i = 0; _complex_specs[i].method != -1; i++)
            index[_complex_specs[i].method] = i;
        initialized = 1;
    }

    return (method >= 0 && method < GR_METHOD_TAB_SIZE && index[method] >= 0) ? _complex_specs + index[method] : NULL;
}

/* f(z) = g(1/z) for large z: acot z = atan(1/z), etc. */
/* inverse functions: exact 1/x feeding the tiny-argument code (e.g. acot(x)
   ~ 1/x for large x), then the generic path */
DECIMAL_DRIVER int
_deccfloat_unary_inv(deccfloat_t res, const deccfloat_t x, int method, int ikind, const ctail_info * tail, gr_ctx_t ctx)
{
    if (!DECFLOAT_IS_ZERO(&x->re) && !DECFLOAT_IS_ZERO(&x->im) && _deccfloat_is_finite(x) && PREC(ctx) != DECIMAL_PREC_EXACT)
    {
        deccfloat_t w;
        int r = 0;
        deccfloat_init(w, ctx);
        if (_deccfloat_inv_exact(w, x, ctx) == GR_SUCCESS)
            r = _deccfloat_tiny_arg(res, w, tail, PREC(ctx), RND(ctx), RND_IM(ctx), ctx);
        deccfloat_clear(w, ctx);
        if (r == 3) return GR_SUCCESS;
        if (r == -1) return GR_UNABLE;
        if (r != 0)
        {
            deccfloat_srcptr a[1];
            deccfloat_ptr rr[1];
            a[0] = x; rr[0] = res;
            return _deccfloat_acb_gr_masked(rr, 1, a, 1, 0, 0, NULL, method, 3 & ~r, ctx);
        }
    }
    return _deccfloat_unary(res, x, method, RK_FN, ikind, NULL, ctx);
}

/* the unary function given by method, componentwise correctly rounded */
int
_deccfloat_unary_method(deccfloat_t res, const deccfloat_t x, int method, gr_ctx_t ctx)
{
    const complex_spec * s = _complex_spec(method);

    if (s == NULL)
        return _deccfloat_unary(res, x, method, RK_NONE, IM_NONE, NULL, ctx);
    if (s->inv)
        return _deccfloat_unary_inv(res, x, method, s->ikind, s->tail, ctx);
    return _deccfloat_unary(res, x, method, s->rkind, s->ikind, s->tail, ctx);
}

int
_deccball_unary_method(deccball_t res, const deccball_t x, int method, gr_ctx_t ctx)
{
    const complex_spec * s = _complex_spec(method);
    return _deccball_unary(res, x, method, s != NULL && s->rkind != RK_NONE && s->rkind != RK_EXP_PI_I, ctx);
}

/* public elementary functions */
#define DEF_CUNARY(name, METHOD) \
int deccfloat_##name(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx) \
{ \
    return _deccfloat_unary_method(res, x, GR_METHOD_##METHOD, ctx); \
} \
int deccball_##name(deccball_t res, const deccball_t x, gr_ctx_t ctx) \
{ \
    return _deccball_unary_method(res, x, GR_METHOD_##METHOD, ctx); \
}

DEF_CUNARY(exp, EXP)
DEF_CUNARY(expm1, EXPM1)
DEF_CUNARY(log, LOG)
DEF_CUNARY(log1p, LOG1P)
DEF_CUNARY(sin, SIN)
DEF_CUNARY(cos, COS)
DEF_CUNARY(tan, TAN)
DEF_CUNARY(asin, ASIN)
DEF_CUNARY(acos, ACOS)
DEF_CUNARY(atan, ATAN)
DEF_CUNARY(sinh, SINH)
DEF_CUNARY(cosh, COSH)
DEF_CUNARY(tanh, TANH)
DEF_CUNARY(asinh, ASINH)
DEF_CUNARY(acosh, ACOSH)
DEF_CUNARY(atanh, ATANH)

/* unary functions with a flag: the real function on the real axis */
DECIMAL_DRIVER int
_deccfloat_unary_flag(deccfloat_t res, const deccfloat_t x, int flag, int method, gr_ctx_t ctx)
{
    deccfloat_srcptr a[1];
    deccfloat_ptr rr[1];
    if (_deccfloat_is_nan(x)) return deccfloat_nan(res, ctx);
    if (DECFLOAT_IS_ZERO(&x->im))
    {
        decfloat_t t;
        int status;
        decfloat_init(t, ctx);
        status = _decfloat_arb_args(t, &x->re, NULL, NULL, NULL, flag, DECIMAL_ARGSPEC(method, 1, 1), ctx);
        if (status == GR_SUCCESS)
        {
            decfloat_swap(&res->re, t, ctx);
            decfloat_zero(&res->im, ctx);
        }
        decfloat_clear(t, ctx);
        if (status == GR_SUCCESS) return status;
    }
    a[0] = x; rr[0] = res;
    return _deccfloat_acb_gr(rr, 1, a, 1, 1, flag, method, ctx);
}

DECIMAL_DRIVER int
_deccball_unary_flag(deccball_t res, const deccball_t x, int flag, int method, gr_ctx_t ctx)
{
    deccball_srcptr a[1];
    deccball_ptr rr[1];
    if (_deccball_is_real(x, ctx))
    {
        decball_t t;
        int status;
        decball_init(t, ctx);
        status = _decball_arb_args(t, &x->re, NULL, NULL, NULL, flag, DECIMAL_ARGSPEC(method, 1, 1), ctx);
        if (status == GR_SUCCESS)
        {
            decball_swap(&res->re, t, ctx);
            decball_zero(&res->im, ctx);
        }
        decball_clear(t, ctx);
        if (status == GR_SUCCESS) return status;
    }
    a[0] = x; rr[0] = res;
    return _deccball_acb_gr(rr, 1, a, 1, 1, flag, method, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Functions of several variables                                         */
/* ------------------------------------------------------------------------- */

/* all arguments real: use the real function; GR_DOMAIN means "use acb" */
DECIMAL_DRIVER int
_deccfloat_multi(deccfloat_ptr * res, slong nres, deccfloat_srcptr * args, slong nargs, int has_flag, int flag, int real, int method, gr_ctx_t ctx)
{
    slong i;
    int all_real = real;

    for (i = 0; i < nargs; i++)
    {
        if (_deccfloat_is_nan(args[i]))
        {
            int status = GR_SUCCESS;
            for (i = 0; i < nres; i++)
                status |= deccfloat_nan(res[i], ctx);
            return status;
        }
        if (!DECFLOAT_IS_ZERO(&args[i]->im))
            all_real = 0;
    }

    if (all_real)
    {
        decfloat_struct t[MAX_RES];
        decfloat_ptr tp[MAX_RES];
        decfloat_srcptr ap[MAX_ARGS];
        int status;

        for (i = 0; i < nres; i++)
        {
            decfloat_init(t + i, ctx);
            tp[i] = t + i;
        }
        for (i = 0; i < nargs; i++)
            ap[i] = &args[i]->re;

        if (nres == 2 && (method == GR_METHOD_SIN_COS || method == GR_METHOD_SIN_COS_PI || method == GR_METHOD_SINH_COSH))
            status = _decfloat_unary2_method(tp[0], tp[1], ap[0],
                (method == GR_METHOD_SIN_COS) ? GR_METHOD_SIN : (method == GR_METHOD_SIN_COS_PI) ? GR_METHOD_SIN_PI : GR_METHOD_SINH,
                (method == GR_METHOD_SIN_COS) ? GR_METHOD_COS : (method == GR_METHOD_SIN_COS_PI) ? GR_METHOD_COS_PI : GR_METHOD_COSH, ctx);
        else
            status = _decfloat_arb_gr(tp, nres, ap, nargs, has_flag, flag, method, ctx);

        if (status == GR_SUCCESS)
        {
            for (i = 0; i < nres; i++)
            {
                decfloat_swap(&res[i]->re, t + i, ctx);
                decfloat_zero(&res[i]->im, ctx);
            }
        }

        for (i = 0; i < nres; i++)
            decfloat_clear(t + i, ctx);

        if (status == GR_SUCCESS)
            return status;
    }

    return _deccfloat_acb_gr(res, nres, args, nargs, has_flag, flag, method, ctx);
}

DECIMAL_DRIVER int
_deccball_multi(deccball_ptr * res, slong nres, deccball_srcptr * args, slong nargs, int has_flag, int flag, int real, int method, gr_ctx_t ctx)
{
    slong i;
    int all_real = real;

    for (i = 0; i < nargs; i++)
        if (!_deccball_is_real(args[i], ctx))
            all_real = 0;

    if (all_real)
    {
        decball_struct t[MAX_RES];
        decball_ptr tp[MAX_RES];
        decball_srcptr ap[MAX_ARGS];
        int status;

        for (i = 0; i < nres; i++)
        {
            decball_init(t + i, ctx);
            tp[i] = t + i;
        }
        for (i = 0; i < nargs; i++)
            ap[i] = &args[i]->re;

        status = _decball_arb_gr(tp, nres, ap, nargs, has_flag, flag, method, ctx);

        if (status == GR_SUCCESS)
        {
            for (i = 0; i < nres; i++)
            {
                decball_swap(&res[i]->re, t + i, ctx);
                decball_zero(&res[i]->im, ctx);
            }
        }

        for (i = 0; i < nres; i++)
            decball_clear(t + i, ctx);

        if (status == GR_SUCCESS)
            return status;
    }

    return _deccball_acb_gr(res, nres, args, nargs, has_flag, flag, method, ctx);
}

/* single-result drivers taking up to four operands directly:
   spec = DECIMAL_ARGSPEC(method, nargs, has_flag) */

DECIMAL_DRIVER int
_deccfloat_multi_args(deccfloat_ptr res, deccfloat_srcptr x, deccfloat_srcptr y, deccfloat_srcptr z, deccfloat_srcptr w, int flag, int real, int spec, gr_ctx_t ctx)
{
    deccfloat_srcptr a[4];
    deccfloat_ptr rr[1];
    a[0] = x; a[1] = y; a[2] = z; a[3] = w; rr[0] = res;
    return _deccfloat_multi(rr, 1, a, (spec / 2) % 8, spec % 2, flag, real, spec / 16, ctx);
}

DECIMAL_DRIVER int
_deccball_multi_args(deccball_ptr res, deccball_srcptr x, deccball_srcptr y, deccball_srcptr z, deccball_srcptr w, int flag, int real, int spec, gr_ctx_t ctx)
{
    deccball_srcptr a[4];
    deccball_ptr rr[1];
    a[0] = x; a[1] = y; a[2] = z; a[3] = w; rr[0] = res;
    return _deccball_multi(rr, 1, a, (spec / 2) % 8, spec % 2, flag, real, spec / 16, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Generic-ring methods for all four decimal types                       */
/* ------------------------------------------------------------------------- */

/*
    The special functions (other than the most common elementary functions,
    which also have public type-specific versions) are only available
    through the generic-ring interface (gr_gamma, gr_bessel_j, ...). Each
    method below handles all four decimal types, according to the context.
*/

#define WHICH(ctx) DECIMAL_CTX_WHICH(ctx)

DECIMAL_DRIVER int
_decimal_unary(gr_ptr res, gr_srcptr x, int method, gr_ctx_t ctx)
{
    switch (WHICH(ctx))
    {
        case DECIMAL_CTX_FLOAT: return _decfloat_unary_method(res, x, method, ctx);
        case DECIMAL_CTX_BALL: return _decball_unary_fn(res, x, method, ctx);
        case DECIMAL_CTX_CFLOAT: return _deccfloat_unary_method(res, x, method, ctx);
        default: return _deccball_unary_method(res, x, method, ctx);
    }
}

/* functions of up to four arguments with an optional flag; real: the
   function is real on the real line (else complex contexts only) */
DECIMAL_DRIVER int
_decimal_multi_args(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_srcptr w, int flag, int real, int spec, gr_ctx_t ctx)
{
    switch (WHICH(ctx))
    {
        case DECIMAL_CTX_FLOAT: return _decfloat_arb_args(res, x, y, z, w, flag, spec, ctx);
        case DECIMAL_CTX_BALL: return _decball_arb_args(res, x, y, z, w, flag, spec, ctx);
        case DECIMAL_CTX_CFLOAT: return _deccfloat_multi_args(res, x, y, z, w, flag, real, spec, ctx);
        default: return _deccball_multi_args(res, x, y, z, w, flag, real, spec, ctx);
    }
}

/* several results */
DECIMAL_DRIVER int
_decimal_multi_res(gr_ptr * res, slong nres, gr_srcptr x, int has_flag, int flag, int method, gr_ctx_t ctx)
{
    gr_srcptr a[1];
    a[0] = x;

    switch (WHICH(ctx))
    {
        case DECIMAL_CTX_FLOAT:
            if (nres == 2 && (method == GR_METHOD_SIN_COS || method == GR_METHOD_SIN_COS_PI || method == GR_METHOD_SINH_COSH))
                return _decfloat_unary2_method(res[0], res[1], x,
                    (method == GR_METHOD_SIN_COS) ? GR_METHOD_SIN : (method == GR_METHOD_SIN_COS_PI) ? GR_METHOD_SIN_PI : GR_METHOD_SINH,
                    (method == GR_METHOD_SIN_COS) ? GR_METHOD_COS : (method == GR_METHOD_SIN_COS_PI) ? GR_METHOD_COS_PI : GR_METHOD_COSH, ctx);
            return _decfloat_arb_gr((decfloat_ptr *) res, nres, (decfloat_srcptr *) a, 1, has_flag, flag, method, ctx);
        case DECIMAL_CTX_BALL:
            return _decball_arb_gr((decball_ptr *) res, nres, (decball_srcptr *) a, 1, has_flag, flag, method, ctx);
        case DECIMAL_CTX_CFLOAT:
            return _deccfloat_multi((deccfloat_ptr *) res, nres, (deccfloat_srcptr *) a, 1, has_flag, flag, 1, method, ctx);
        default:
            return _deccball_multi((deccball_ptr *) res, nres, (deccball_srcptr *) a, 1, has_flag, flag, 1, method, ctx);
    }
}

/* constants (real) */
DECIMAL_DRIVER int
_decimal_constant(gr_ptr res, int method, gr_ctx_t ctx)
{
    switch (WHICH(ctx))
    {
        case DECIMAL_CTX_FLOAT: return _decfloat_constant(res, method, ctx);
        case DECIMAL_CTX_BALL: return _decball_constant(res, method, ctx);
        case DECIMAL_CTX_CFLOAT:
            decfloat_zero(DECCFLOAT_IMAGREF((deccfloat_ptr) res), ctx);
            return _decfloat_constant(DECCFLOAT_REALREF((deccfloat_ptr) res), method, ctx);
        default:
            decball_zero(DECCBALL_IMAGREF((deccball_ptr) res), ctx);
            return _decball_constant(DECCBALL_REALREF((deccball_ptr) res), method, ctx);
    }
}

/* f(0) = 0 exactly for the Fresnel integrals (arb may not give an exact zero) */
static int
_decimal_zero_at_zero(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    switch (WHICH(ctx))
    {
        case DECIMAL_CTX_FLOAT:
            return DECFLOAT_IS_ZERO((decfloat_srcptr) x) && decfloat_zero(res, ctx) == GR_SUCCESS;
        case DECIMAL_CTX_CFLOAT:
            return DECFLOAT_IS_ZERO(DECCFLOAT_REALREF((deccfloat_srcptr) x)) && DECFLOAT_IS_ZERO(DECCFLOAT_IMAGREF((deccfloat_srcptr) x)) && deccfloat_zero(res, ctx) == GR_SUCCESS;
        default:
            return 0;
    }
}

#define DEF_GR_UNARY(name, METHOD) \
static int _decimal_##name(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) \
{ return _decimal_unary(res, x, GR_METHOD_##METHOD, ctx); }

#define DEF_GR_UNARY_FLAG(name, METHOD, ZERO_AT_ZERO) \
static int _decimal_##name(gr_ptr res, gr_srcptr x, int flag, gr_ctx_t ctx) \
{ \
    if (ZERO_AT_ZERO && _decimal_zero_at_zero(res, x, ctx)) return GR_SUCCESS; \
    switch (WHICH(ctx)) \
    { \
        case DECIMAL_CTX_CFLOAT: return _deccfloat_unary_flag(res, x, flag, GR_METHOD_##METHOD, ctx); \
        case DECIMAL_CTX_CBALL: return _deccball_unary_flag(res, x, flag, GR_METHOD_##METHOD, ctx); \
        default: return _decimal_multi_args(res, x, NULL, NULL, NULL, flag, 1, DECIMAL_ARGSPEC(GR_METHOD_##METHOD, 1, 1), ctx); \
    } \
}

#define DEF_GR_BINARY(name, METHOD, REAL) \
static int _decimal_##name(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) \
{ return _decimal_multi_args(res, x, y, NULL, NULL, 0, REAL, DECIMAL_ARGSPEC(GR_METHOD_##METHOD, 2, 0), ctx); }

#define DEF_GR_BINARY_FLAG(name, METHOD, REAL) \
static int _decimal_##name(gr_ptr res, gr_srcptr x, gr_srcptr y, int flag, gr_ctx_t ctx) \
{ return _decimal_multi_args(res, x, y, NULL, NULL, flag, REAL, DECIMAL_ARGSPEC(GR_METHOD_##METHOD, 2, 1), ctx); }

#define DEF_GR_TERNARY(name, METHOD, REAL) \
static int _decimal_##name(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) \
{ return _decimal_multi_args(res, x, y, z, NULL, 0, REAL, DECIMAL_ARGSPEC(GR_METHOD_##METHOD, 3, 0), ctx); }

#define DEF_GR_TERNARY_FLAG(name, METHOD, REAL) \
static int _decimal_##name(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int flag, gr_ctx_t ctx) \
{ return _decimal_multi_args(res, x, y, z, NULL, flag, REAL, DECIMAL_ARGSPEC(GR_METHOD_##METHOD, 3, 1), ctx); }

#define DEF_GR_CONSTANT(name, METHOD) \
static int _decimal_##name(gr_ptr res, gr_ctx_t ctx) \
{ return _decimal_constant(res, GR_METHOD_##METHOD, ctx); }

DEF_GR_UNARY(exp, EXP)
DEF_GR_UNARY(expm1, EXPM1)
DEF_GR_UNARY(exp2, EXP2)
DEF_GR_UNARY(exp10, EXP10)
DEF_GR_UNARY(log, LOG)
DEF_GR_UNARY(log1p, LOG1P)
DEF_GR_UNARY(log2, LOG2)
DEF_GR_UNARY(log10, LOG10)
DEF_GR_UNARY(sin, SIN)
DEF_GR_UNARY(cos, COS)
DEF_GR_UNARY(tan, TAN)
DEF_GR_UNARY(cot, COT)
DEF_GR_UNARY(sec, SEC)
DEF_GR_UNARY(csc, CSC)
DEF_GR_UNARY(sinc, SINC)
DEF_GR_UNARY(sin_pi, SIN_PI)
DEF_GR_UNARY(cos_pi, COS_PI)
DEF_GR_UNARY(tan_pi, TAN_PI)
DEF_GR_UNARY(cot_pi, COT_PI)
DEF_GR_UNARY(sec_pi, SEC_PI)
DEF_GR_UNARY(csc_pi, CSC_PI)
DEF_GR_UNARY(sinc_pi, SINC_PI)
DEF_GR_UNARY(asin, ASIN)
DEF_GR_UNARY(acos, ACOS)
DEF_GR_UNARY(atan, ATAN)
DEF_GR_UNARY(acot, ACOT)
DEF_GR_UNARY(asec, ASEC)
DEF_GR_UNARY(acsc, ACSC)
DEF_GR_UNARY(asin_pi, ASIN_PI)
DEF_GR_UNARY(acos_pi, ACOS_PI)
DEF_GR_UNARY(atan_pi, ATAN_PI)
DEF_GR_UNARY(acot_pi, ACOT_PI)
DEF_GR_UNARY(asec_pi, ASEC_PI)
DEF_GR_UNARY(acsc_pi, ACSC_PI)
DEF_GR_UNARY(sinh, SINH)
DEF_GR_UNARY(cosh, COSH)
DEF_GR_UNARY(tanh, TANH)
DEF_GR_UNARY(coth, COTH)
DEF_GR_UNARY(sech, SECH)
DEF_GR_UNARY(csch, CSCH)
DEF_GR_UNARY(asinh, ASINH)
DEF_GR_UNARY(acosh, ACOSH)
DEF_GR_UNARY(atanh, ATANH)
DEF_GR_UNARY(acoth, ACOTH)
DEF_GR_UNARY(asech, ASECH)
DEF_GR_UNARY(acsch, ACSCH)
DEF_GR_UNARY(lambertw, LAMBERTW)
DEF_GR_UNARY(gamma, GAMMA)
DEF_GR_UNARY(rgamma, RGAMMA)
DEF_GR_UNARY(lgamma, LGAMMA)
DEF_GR_UNARY(digamma, DIGAMMA)
DEF_GR_UNARY(barnes_g, BARNES_G)
DEF_GR_UNARY(log_barnes_g, LOG_BARNES_G)
DEF_GR_UNARY(zeta, ZETA)
DEF_GR_UNARY(erf, ERF)
DEF_GR_UNARY(erfc, ERFC)
DEF_GR_UNARY(erfi, ERFI)
DEF_GR_UNARY(erfinv, ERFINV)
DEF_GR_UNARY(erfcinv, ERFCINV)
DEF_GR_UNARY(exp_integral_ei, EXP_INTEGRAL_EI)
DEF_GR_UNARY(sin_integral, SIN_INTEGRAL)
DEF_GR_UNARY(cos_integral, COS_INTEGRAL)
DEF_GR_UNARY(sinh_integral, SINH_INTEGRAL)
DEF_GR_UNARY(cosh_integral, COSH_INTEGRAL)
DEF_GR_UNARY(dilog, DILOG)
DEF_GR_UNARY(agm1, AGM1)
DEF_GR_UNARY(airy_ai, AIRY_AI)
DEF_GR_UNARY(airy_bi, AIRY_BI)
DEF_GR_UNARY(airy_ai_prime, AIRY_AI_PRIME)
DEF_GR_UNARY(airy_bi_prime, AIRY_BI_PRIME)
DEF_GR_UNARY(exp_pi_i, EXP_PI_I)
DEF_GR_UNARY(log_pi_i, LOG_PI_I)
DEF_GR_UNARY(dirichlet_eta, DIRICHLET_ETA)
DEF_GR_UNARY(riemann_xi, RIEMANN_XI)
DEF_GR_UNARY(modular_j, MODULAR_J)
DEF_GR_UNARY(modular_lambda, MODULAR_LAMBDA)
DEF_GR_UNARY(modular_delta, MODULAR_DELTA)
DEF_GR_UNARY(dedekind_eta, DEDEKIND_ETA)
DEF_GR_UNARY(elliptic_k, ELLIPTIC_K)
DEF_GR_UNARY(elliptic_e, ELLIPTIC_E)

DEF_GR_UNARY_FLAG(fresnel_s, FRESNEL_S, 1)
DEF_GR_UNARY_FLAG(fresnel_c, FRESNEL_C, 1)
DEF_GR_UNARY_FLAG(log_integral, LOG_INTEGRAL, 0)

DEF_GR_BINARY(agm, AGM, 1)
DEF_GR_BINARY(rising, RISING, 1)
DEF_GR_BINARY(bessel_j, BESSEL_J, 1)
DEF_GR_BINARY(bessel_y, BESSEL_Y, 1)
DEF_GR_BINARY(bessel_i, BESSEL_I, 1)
DEF_GR_BINARY(bessel_k, BESSEL_K, 1)
DEF_GR_BINARY(bessel_i_scaled, BESSEL_I_SCALED, 1)
DEF_GR_BINARY(bessel_k_scaled, BESSEL_K_SCALED, 1)
DEF_GR_BINARY(polylog, POLYLOG, 1)
DEF_GR_BINARY(hurwitz_zeta, HURWITZ_ZETA, 1)
DEF_GR_BINARY(exp_integral, EXP_INTEGRAL, 1)
DEF_GR_BINARY(chebyshev_t, CHEBYSHEV_T, 1)
DEF_GR_BINARY(chebyshev_u, CHEBYSHEV_U, 1)
DEF_GR_BINARY(hermite_h, HERMITE_H, 1)
DEF_GR_BINARY(polygamma, POLYGAMMA, 0)
DEF_GR_BINARY(elliptic_pi, ELLIPTIC_PI, 0)

DEF_GR_BINARY_FLAG(gamma_upper, GAMMA_UPPER, 1)
DEF_GR_BINARY_FLAG(gamma_lower, GAMMA_LOWER, 1)
DEF_GR_BINARY_FLAG(hypgeom_0f1, HYPGEOM_0F1, 1)
DEF_GR_BINARY_FLAG(elliptic_f, ELLIPTIC_F, 0)
DEF_GR_BINARY_FLAG(elliptic_e_inc, ELLIPTIC_E_INC, 0)

DEF_GR_TERNARY(coulomb_f, COULOMB_F, 1)
DEF_GR_TERNARY(coulomb_g, COULOMB_G, 1)
DEF_GR_TERNARY(gegenbauer_c, GEGENBAUER_C, 1)
DEF_GR_TERNARY(laguerre_l, LAGUERRE_L, 1)
DEF_GR_TERNARY(lerch_phi, LERCH_PHI, 0)

DEF_GR_TERNARY_FLAG(beta_lower, BETA_LOWER, 1)
DEF_GR_TERNARY_FLAG(hypgeom_1f1, HYPGEOM_1F1, 1)
DEF_GR_TERNARY_FLAG(hypgeom_u, HYPGEOM_U, 1)
DEF_GR_TERNARY_FLAG(legendre_p, LEGENDRE_P, 1)
DEF_GR_TERNARY_FLAG(legendre_q, LEGENDRE_Q, 1)

static int _decimal_hypgeom_2f1(gr_ptr res, gr_srcptr a, gr_srcptr b, gr_srcptr c, gr_srcptr z, int flag, gr_ctx_t ctx)
{ return _decimal_multi_args(res, a, b, c, z, flag, 1, DECIMAL_ARGSPEC(GR_METHOD_HYPGEOM_2F1, 4, 1), ctx); }
static int _decimal_jacobi_p(gr_ptr res, gr_srcptr n, gr_srcptr a, gr_srcptr b, gr_srcptr z, gr_ctx_t ctx)
{ return _decimal_multi_args(res, n, a, b, z, 0, 1, DECIMAL_ARGSPEC(GR_METHOD_JACOBI_P, 4, 0), ctx); }

static int _decimal_sin_cos(gr_ptr r1, gr_ptr r2, gr_srcptr x, gr_ctx_t ctx)
{ gr_ptr r[2]; r[0] = r1; r[1] = r2; return _decimal_multi_res(r, 2, x, 0, 0, GR_METHOD_SIN_COS, ctx); }
static int _decimal_sin_cos_pi(gr_ptr r1, gr_ptr r2, gr_srcptr x, gr_ctx_t ctx)
{ gr_ptr r[2]; r[0] = r1; r[1] = r2; return _decimal_multi_res(r, 2, x, 0, 0, GR_METHOD_SIN_COS_PI, ctx); }
static int _decimal_sinh_cosh(gr_ptr r1, gr_ptr r2, gr_srcptr x, gr_ctx_t ctx)
{ gr_ptr r[2]; r[0] = r1; r[1] = r2; return _decimal_multi_res(r, 2, x, 0, 0, GR_METHOD_SINH_COSH, ctx); }
static int _decimal_fresnel(gr_ptr r1, gr_ptr r2, gr_srcptr x, int normalized, gr_ctx_t ctx)
{ gr_ptr r[2]; r[0] = r1; r[1] = r2; return _decimal_multi_res(r, 2, x, 1, normalized, GR_METHOD_FRESNEL, ctx); }
static int _decimal_airy(gr_ptr ai, gr_ptr aip, gr_ptr bi, gr_ptr bip, gr_srcptr x, gr_ctx_t ctx)
{ gr_ptr r[4]; r[0] = ai; r[1] = aip; r[2] = bi; r[3] = bip; return _decimal_multi_res(r, 4, x, 0, 0, GR_METHOD_AIRY, ctx); }

DEF_GR_CONSTANT(pi, PI)
DEF_GR_CONSTANT(euler, EULER)
DEF_GR_CONSTANT(catalan, CATALAN)
DEF_GR_CONSTANT(khinchin, KHINCHIN)
DEF_GR_CONSTANT(glaisher, GLAISHER)

/* functions of integers and rationals (real) */
#define DEF_GR_REAL_OF(name, T) \
static int _decimal_##name(gr_ptr res, T n, gr_ctx_t ctx) \
{ \
    switch (WHICH(ctx)) \
    { \
        case DECIMAL_CTX_FLOAT: return _decfloat_##name(res, n, ctx); \
        case DECIMAL_CTX_BALL: return _decball_##name(res, n, ctx); \
        case DECIMAL_CTX_CFLOAT: \
            decfloat_zero(DECCFLOAT_IMAGREF((deccfloat_ptr) res), ctx); \
            return _decfloat_##name(DECCFLOAT_REALREF((deccfloat_ptr) res), n, ctx); \
        default: \
            decball_zero(DECCBALL_IMAGREF((deccball_ptr) res), ctx); \
            return _decball_##name(DECCBALL_REALREF((deccball_ptr) res), n, ctx); \
    } \
}

DEF_GR_REAL_OF(gamma_fmpz, const fmpz_t)
DEF_GR_REAL_OF(gamma_fmpq, const fmpq_t)
DEF_GR_REAL_OF(fac_ui, ulong)
DEF_GR_REAL_OF(fac_fmpz, const fmpz_t)

/* rising factorial with an integer argument: exact when feasible */
static int
_deccfloat_rising_ui(deccfloat_t res, const deccfloat_t x, ulong n, gr_ctx_t ctx)
{
    extra_arg extra = {1, n, NULL, NULL};
    deccfloat_srcptr a[1];
    deccfloat_ptr rr[1];

    if (DECFLOAT_IS_ZERO(&x->im))
    {
        decfloat_t t;
        int status;
        decfloat_init(t, ctx);
        status = _decfloat_rising_ui(t, &x->re, n, ctx);
        if (status == GR_SUCCESS)
        {
            decfloat_swap(&res->re, t, ctx);
            decfloat_zero(&res->im, ctx);
        }
        decfloat_clear(t, ctx);
        return status;
    }

    if (n == 0)
        return deccfloat_one(res, ctx);

    if (_deccfloat_is_finite(x) && n <= 1000 &&
        (decfloat_digits(&x->re, ctx) + decfloat_digits(&x->im, ctx) + 4
            + (DECFLOAT_IS_ZERO(&x->re) ? 0 : FLINT_MIN(FLINT_ABS(_decfloat_val10_clamped(&x->re, ctx)), 1000000))
            + (DECFLOAT_IS_ZERO(&x->im) ? 0 : FLINT_MIN(FLINT_ABS(_decfloat_val10_clamped(&x->im, ctx)), 1000000))) * n <= 100000)
    {
        deccfloat_t t, u;
        ulong k;
        int status = GR_SUCCESS;

        deccfloat_init(t, ctx);
        deccfloat_init(u, ctx);
        deccfloat_one(t, ctx);

        for (k = 0; k < n && status == GR_SUCCESS; k++)
        {
            status = decfloat_set_round_ui(&u->re, k, DECIMAL_PREC_EXACT, EXACT_RND, ctx);
            decfloat_zero(&u->im, ctx);
            status |= _deccfloat_add_exact(u, x, u, 0, ctx);
            if (status == GR_SUCCESS)
                status = _deccfloat_mul_exact(t, t, u, ctx);
        }

        if (status == GR_SUCCESS)
            status = deccfloat_set(res, t, ctx);

        deccfloat_clear(t, ctx);
        deccfloat_clear(u, ctx);
        return status;
    }

    a[0] = x;
    rr[0] = res;
    return _deccfloat_acb_gr_extra(rr, 1, a, 1, 0, 0, &extra, GR_METHOD_RISING_UI, ctx);
}

static int
_deccball_rising_ui(deccball_t res, const deccball_t x, ulong n, gr_ctx_t ctx)
{
    extra_arg extra = {1, n, NULL, NULL};
    deccball_srcptr a[1];
    deccball_ptr rr[1];

    if (_deccball_is_real(x, ctx))
    {
        decball_zero(&res->im, ctx);
        return _decball_rising_ui(&res->re, &x->re, n, ctx);
    }

    a[0] = x;
    rr[0] = res;
    return _deccball_acb_gr_extra(rr, 1, a, 1, 0, 0, &extra, GR_METHOD_RISING_UI, ctx);
}

static int
_deccfloat_lambertw_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t k, gr_ctx_t ctx)
{
    extra_arg extra = {2, 0, k, NULL};
    deccfloat_srcptr a[1];
    deccfloat_ptr rr[1];

    if (fmpz_is_zero(k))
        return _deccfloat_unary_method(res, x, GR_METHOD_LAMBERTW, ctx);

    a[0] = x;
    rr[0] = res;
    return _deccfloat_acb_gr_extra(rr, 1, a, 1, 0, 0, &extra, GR_METHOD_LAMBERTW_FMPZ, ctx);
}

static int
_deccball_lambertw_fmpz(deccball_t res, const deccball_t x, const fmpz_t k, gr_ctx_t ctx)
{
    extra_arg extra = {2, 0, k, NULL};
    deccball_srcptr a[1];
    deccball_ptr rr[1];

    if (fmpz_is_zero(k))
        return _deccball_unary_method(res, x, GR_METHOD_LAMBERTW, ctx);

    a[0] = x;
    rr[0] = res;
    return _deccball_acb_gr_extra(rr, 1, a, 1, 0, 0, &extra, GR_METHOD_LAMBERTW_FMPZ, ctx);
}


#define DEF_GR_WITH_INT(name, T) \
static int _decimal_##name(gr_ptr res, gr_srcptr x, T n, gr_ctx_t ctx) \
{ \
    switch (WHICH(ctx)) \
    { \
        case DECIMAL_CTX_FLOAT: return _decfloat_##name(res, x, n, ctx); \
        case DECIMAL_CTX_BALL: return _decball_##name(res, x, n, ctx); \
        case DECIMAL_CTX_CFLOAT: return _deccfloat_##name(res, x, n, ctx); \
        default: return _deccball_##name(res, x, n, ctx); \
    } \
}

DEF_GR_WITH_INT(rising_ui, ulong)
DEF_GR_WITH_INT(lambertw_fmpz, const fmpz_t)

/* method tables shared by the four decimal types (applied on top of the
   type-specific tables): methods of all types, and of the complex types */
#define M(METHOD, name) {GR_METHOD_##METHOD, (gr_funcptr) _decimal_##name}

gr_method_tab_input _decimal_special_methods_input[] =
{
    M(EXP, exp), M(EXPM1, expm1), M(EXP2, exp2), M(EXP10, exp10),
    M(LOG, log), M(LOG1P, log1p), M(LOG2, log2), M(LOG10, log10),
    M(SIN, sin), M(COS, cos), M(TAN, tan), M(COT, cot), M(SEC, sec), M(CSC, csc), M(SINC, sinc),
    M(SIN_PI, sin_pi), M(COS_PI, cos_pi), M(TAN_PI, tan_pi), M(COT_PI, cot_pi), M(SEC_PI, sec_pi),
    M(CSC_PI, csc_pi), M(SINC_PI, sinc_pi),
    M(ASIN, asin), M(ACOS, acos), M(ATAN, atan), M(ACOT, acot), M(ASEC, asec), M(ACSC, acsc),
    M(ASIN_PI, asin_pi), M(ACOS_PI, acos_pi), M(ATAN_PI, atan_pi), M(ACOT_PI, acot_pi),
    M(ASEC_PI, asec_pi), M(ACSC_PI, acsc_pi),
    M(SINH, sinh), M(COSH, cosh), M(TANH, tanh), M(COTH, coth), M(SECH, sech), M(CSCH, csch),
    M(ASINH, asinh), M(ACOSH, acosh), M(ATANH, atanh), M(ACOTH, acoth), M(ASECH, asech), M(ACSCH, acsch),
    M(SIN_COS, sin_cos), M(SIN_COS_PI, sin_cos_pi), M(SINH_COSH, sinh_cosh),
    M(LAMBERTW, lambertw), M(LAMBERTW_FMPZ, lambertw_fmpz),
    M(GAMMA, gamma), M(RGAMMA, rgamma), M(LGAMMA, lgamma), M(DIGAMMA, digamma),
    M(GAMMA_FMPZ, gamma_fmpz), M(GAMMA_FMPQ, gamma_fmpq), M(FAC_UI, fac_ui), M(FAC_FMPZ, fac_fmpz),
    M(RISING_UI, rising_ui), M(RISING, rising),
    M(BARNES_G, barnes_g), M(LOG_BARNES_G, log_barnes_g), M(ZETA, zeta), M(HURWITZ_ZETA, hurwitz_zeta),
    M(POLYLOG, polylog),
    M(ERF, erf), M(ERFC, erfc), M(ERFI, erfi), M(ERFINV, erfinv), M(ERFCINV, erfcinv),
    M(FRESNEL, fresnel), M(FRESNEL_S, fresnel_s), M(FRESNEL_C, fresnel_c),
    M(EXP_INTEGRAL, exp_integral), M(EXP_INTEGRAL_EI, exp_integral_ei),
    M(SIN_INTEGRAL, sin_integral), M(COS_INTEGRAL, cos_integral),
    M(SINH_INTEGRAL, sinh_integral), M(COSH_INTEGRAL, cosh_integral), M(LOG_INTEGRAL, log_integral),
    M(DILOG, dilog), M(AGM, agm), M(AGM1, agm1),
    M(AIRY, airy), M(AIRY_AI, airy_ai), M(AIRY_BI, airy_bi), M(AIRY_AI_PRIME, airy_ai_prime),
    M(AIRY_BI_PRIME, airy_bi_prime),
    M(BESSEL_J, bessel_j), M(BESSEL_Y, bessel_y), M(BESSEL_I, bessel_i), M(BESSEL_K, bessel_k),
    M(BESSEL_I_SCALED, bessel_i_scaled), M(BESSEL_K_SCALED, bessel_k_scaled),
    M(GAMMA_UPPER, gamma_upper), M(GAMMA_LOWER, gamma_lower), M(BETA_LOWER, beta_lower),
    M(HYPGEOM_0F1, hypgeom_0f1), M(HYPGEOM_1F1, hypgeom_1f1), M(HYPGEOM_U, hypgeom_u),
    M(HYPGEOM_2F1, hypgeom_2f1),
    M(CHEBYSHEV_T, chebyshev_t), M(CHEBYSHEV_U, chebyshev_u), M(JACOBI_P, jacobi_p),
    M(GEGENBAUER_C, gegenbauer_c), M(LAGUERRE_L, laguerre_l), M(HERMITE_H, hermite_h),
    M(LEGENDRE_P, legendre_p), M(LEGENDRE_Q, legendre_q),
    M(COULOMB_F, coulomb_f), M(COULOMB_G, coulomb_g),
    M(PI, pi), M(EULER, euler), M(CATALAN, catalan), M(KHINCHIN, khinchin), M(GLAISHER, glaisher),
    {0, (gr_funcptr) NULL}
};

gr_method_tab_input _decimal_complex_special_methods_input[] =
{
    M(EXP_PI_I, exp_pi_i), M(LOG_PI_I, log_pi_i),
    M(DIRICHLET_ETA, dirichlet_eta), M(RIEMANN_XI, riemann_xi), M(POLYGAMMA, polygamma),
    M(LERCH_PHI, lerch_phi),
    M(MODULAR_J, modular_j), M(MODULAR_LAMBDA, modular_lambda), M(MODULAR_DELTA, modular_delta),
    M(DEDEKIND_ETA, dedekind_eta),
    M(ELLIPTIC_K, elliptic_k), M(ELLIPTIC_E, elliptic_e), M(ELLIPTIC_PI, elliptic_pi),
    M(ELLIPTIC_F, elliptic_f), M(ELLIPTIC_E_INC, elliptic_e_inc),
    {0, (gr_funcptr) NULL}
};

#undef M

/* public versions of the most common elementary functions */
int deccfloat_sin_cos(deccfloat_t res1, deccfloat_t res2, const deccfloat_t x, gr_ctx_t ctx) { return _decimal_sin_cos(res1, res2, x, ctx); }
int deccball_sin_cos(deccball_t res1, deccball_t res2, const deccball_t x, gr_ctx_t ctx) { return _decimal_sin_cos(res1, res2, x, ctx); }
int deccfloat_sinh_cosh(deccfloat_t res1, deccfloat_t res2, const deccfloat_t x, gr_ctx_t ctx) { return _decimal_sinh_cosh(res1, res2, x, ctx); }
int deccball_sinh_cosh(deccball_t res1, deccball_t res2, const deccball_t x, gr_ctx_t ctx) { return _decimal_sinh_cosh(res1, res2, x, ctx); }
int deccfloat_pi(deccfloat_t res, gr_ctx_t ctx) { return _decimal_pi(res, ctx); }
int deccball_pi(deccball_t res, gr_ctx_t ctx) { return _decimal_pi(res, ctx); }

POP_OPTIONS
