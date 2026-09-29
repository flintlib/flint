/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Elementary and special functions for decimal floats and balls,
    computed through the generic-ring methods of arb.

    Balls: each argument is converted to an arb ball, the arb method is
    called at the working precision, and the results are converted back.

    Floats: the same, inside Ziv's loop, until the rounding of every
    result is decided. Before that, exact special values (such as
    Gamma(5) = 24, sin(pi/6) = 1/2, log10(1000) = 3) and tiny arguments
    (where the exact result lies within a sticky bit of a Taylor partial
    sum) are handled directly, since Ziv's loop cannot terminate on them
    in the directed rounding modes.
*/

#include <math.h>
#include "decimal.h"
#include "mag.h"
#include "fmpq.h"
#include "arith.h"
#include "arf.h"
#include "arb.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_special.h"

/* Not performance-critical: optimize for size. */
PUSH_OPTIONS
OPTIMIZE_OSIZE

/* ------------------------------------------------------------------------- */
/*    Calling arb methods                                                    */
/* ------------------------------------------------------------------------- */

/* an extra integer or rational argument passed straight through */
typedef struct
{
    int type;       /* 0: none, 1: ulong, 2: fmpz, 3: fmpq */
    ulong u;
    const fmpz * z;
    const fmpq * q;
}
extra_arg;

DECIMAL_DRIVER int
_arb_dispatch(arb_ptr * res, slong nres, arb_srcptr * args, slong nargs, int has_flag, int flag, const extra_arg * extra, gr_funcptr fn, gr_ctx_t actx)
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
                flint_throw(FLINT_ERROR, "_arb_dispatch: unsupported signature\n");
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
            flint_throw(FLINT_ERROR, "_arb_dispatch: unsupported signature\n");
    }
}

#define MAX_ARGS 4
#define MAX_RES 4

/* arb method on balls */
DECIMAL_DRIVER int
_decball_arb_gr_extra(decball_ptr * res, slong nres, decball_srcptr * args, slong nargs, int has_flag, int flag, const extra_arg * extra, int method, gr_ctx_t ctx)
{
    gr_ctx_t actx;
    arb_struct a[MAX_ARGS], r[MAX_RES];
    arb_ptr ap[MAX_ARGS], rp[MAX_RES];
    arb_srcptr asp[MAX_ARGS];
    decfloat_srcptr mids[MAX_ARGS];
    slong i, wp;
    int status = GR_SUCCESS;

    for (i = 0; i < nargs; i++)
        mids[i] = &args[i]->mid;
    if (_decimal_arb_args_infeasible(method, mids, NULL, nargs, DECIMAL_CTX_PREC(ctx), ctx))
        return GR_UNABLE;

    wp = _decimal_digits_to_bits(DECIMAL_CTX_PREC(ctx));
    gr_ctx_init_real_arb(actx, wp);

    for (i = 0; i < nargs; i++)
    {
        arb_init(a + i);
        ap[i] = a + i;
        asp[i] = a + i;
        status |= decball_get_arb(a + i, args[i], wp, ctx);
    }
    for (i = 0; i < nres; i++)
    {
        arb_init(r + i);
        rp[i] = r + i;
    }

    if (status == GR_SUCCESS)
        status = _arb_dispatch(rp, nres, asp, nargs, has_flag, flag, extra, actx->methods[method], actx);

    if (status == GR_SUCCESS)
        for (i = 0; i < nres; i++)
            status |= decball_set_arb(res[i], r + i, ctx);

    for (i = 0; i < nargs; i++)
        arb_clear(a + i);
    for (i = 0; i < nres; i++)
        arb_clear(r + i);
    gr_ctx_clear(actx);
    (void) ap;
    return status;
}

/* arb method on floats, correctly rounded with Ziv's strategy */
DECIMAL_DRIVER int
_decfloat_arb_gr_extra(decfloat_ptr * res, slong nres, decfloat_srcptr * args, slong nargs, int has_flag, int flag, const extra_arg * extra, int method, gr_ctx_t ctx)
{
    slong prec = DECIMAL_CTX_PREC(ctx);
    int rnd = DECIMAL_CTX_RND(ctx);
    gr_ctx_t actx, bctx;
    arb_struct a[MAX_ARGS], r[MAX_RES];
    arb_ptr rp[MAX_RES];
    arb_srcptr asp[MAX_ARGS];
    decball_t Y;
    decfloat_struct tmp[MAX_RES];
    slong i, wp, wpbits, wp_max;
    int status = GR_UNABLE, rr;

    for (i = 0; i < nargs; i++)
    {
        if (DECFLOAT_IS_NAN(args[i]))
        {
            for (i = 0; i < nres; i++)
                status = decfloat_nan(res[i], ctx);
            return status;
        }
        DECFLOAT_CHECK_OPERAND(args[i], ctx);
    }

    if (prec == DECIMAL_PREC_EXACT)
        return GR_UNABLE;

    if (_decimal_arb_args_infeasible(method, args, NULL, nargs, prec, ctx))
        return GR_UNABLE;

    _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), prec, DECIMAL_RND_DOWN, 0);
    decimal_ctx_set_rad_prec(bctx, DECMAG_MAX_PREC);
    decball_init(Y, bctx);

    for (i = 0; i < nargs; i++)
    {
        arb_init(a + i);
        asp[i] = a + i;
    }
    for (i = 0; i < nres; i++)
    {
        arb_init(r + i);
        rp[i] = r + i;
        decfloat_init(tmp + i, ctx);
    }

    wp_max = 0;
    if (_decimal_method_is_elementary(method))
        for (i = 0; i < nargs; i++)
            wp_max = FLINT_MAX(wp_max, decfloat_digits(args[i], ctx));
    wp_max = _decimal_ziv_wp_max(prec, wp_max);

    for (wp = prec + 10; wp <= wp_max; wp *= 2)
    {
        wpbits = _decimal_digits_to_bits(wp);
        decimal_ctx_set_prec(bctx, wp);
        gr_ctx_init_real_arb(actx, wpbits);

        status = GR_SUCCESS;
        for (i = 0; i < nargs; i++)
            status |= decfloat_get_arb(a + i, args[i], wpbits, ctx);

        if (status == GR_SUCCESS)
            status = _arb_dispatch(rp, nres, asp, nargs, has_flag, flag, extra, actx->methods[method], actx);

        gr_ctx_clear(actx);

        if (status != GR_SUCCESS)
        {
            /* GR_UNABLE at a low precision can just mean that arb needs
               more precision (for example to rule out a pole) */
            if (status & GR_DOMAIN)
                break;
            continue;
        }

        rr = 1;
        for (i = 0; i < nres && rr == 1; i++)
        {
            if (!arb_is_finite(r + i))
            {
                if (arf_is_pos_inf(arb_midref(r + i)) && mag_is_finite(arb_radref(r + i)))
                    rr = (decfloat_pos_inf(tmp + i, ctx) == GR_SUCCESS) ? 1 : -1;
                else if (arf_is_neg_inf(arb_midref(r + i)) && mag_is_finite(arb_radref(r + i)))
                    rr = (decfloat_neg_inf(tmp + i, ctx) == GR_SUCCESS) ? 1 : -1;
                else
                    rr = 0;   /* intermediate overflow: try again */
                continue;
            }

            if (decball_set_arb(Y, r + i, bctx) != GR_SUCCESS)
            {
                rr = -1;
                break;
            }

            rr = _decfloat_round_ball(tmp + i, Y, prec, rnd, bctx, ctx);
            if (rr == 1)
                rr = (_decfloat_finalize(tmp + i, ctx) == GR_SUCCESS) ? 1 : -1;
        }

        if (rr == 1)
        {
            for (i = 0; i < nres; i++)
                decfloat_swap(res[i], tmp + i, ctx);
            status = GR_SUCCESS;
            break;
        }

        status = GR_UNABLE;
        if (rr == -1)
            break;
    }

    for (i = 0; i < nargs; i++)
        arb_clear(a + i);
    for (i = 0; i < nres; i++)
    {
        arb_clear(r + i);
        decfloat_clear(tmp + i, ctx);
    }
    decball_clear(Y, bctx);
    gr_ctx_clear(bctx);
    return status;
}

int
_decball_arb_gr(decball_ptr * res, slong nres, decball_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx)
{
    return _decball_arb_gr_extra(res, nres, args, nargs, has_flag, flag, NULL, method, ctx);
}

int
_decfloat_arb_gr(decfloat_ptr * res, slong nres, decfloat_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx)
{
    return _decfloat_arb_gr_extra(res, nres, args, nargs, has_flag, flag, NULL, method, ctx);
}

/* Single-result drivers taking up to four operands directly, so that the
   wrappers are plain tail calls: spec = DECIMAL_ARGSPEC(method, nargs, has_flag). */

int
_decfloat_arb_args(decfloat_ptr res, decfloat_srcptr x, decfloat_srcptr y, decfloat_srcptr z, decfloat_srcptr w, int flag, int spec, gr_ctx_t ctx)
{
    decfloat_srcptr a[4];
    decfloat_ptr rr[1];
    a[0] = x; a[1] = y; a[2] = z; a[3] = w; rr[0] = res;
    return _decfloat_arb_gr(rr, 1, a, (spec / 2) % 8, spec % 2, flag, spec / 16, ctx);
}

int
_decball_arb_args(decball_ptr res, decball_srcptr x, decball_srcptr y, decball_srcptr z, decball_srcptr w, int flag, int spec, gr_ctx_t ctx)
{
    decball_srcptr a[4];
    decball_ptr rr[1];
    a[0] = x; a[1] = y; a[2] = z; a[3] = w; rr[0] = res;
    return _decball_arb_gr(rr, 1, a, (spec / 2) % 8, spec % 2, flag, spec / 16, ctx);
}

/*
    Arguments for which an evaluation through arb/acb is not attempted
    (GR_UNABLE is returned), since the arb algorithms would need time or
    memory far beyond what the result is worth:

    * for non-elementary functions, a part of huge magnitude (scientific
      exponent > DECIMAL_ARB_EXP_LIMIT), and for the Barnes G-function also
      of tiny magnitude: the arb algorithms then work with numbers whose
      size is proportional to the exponent;
    * for the polylogarithm, the Hurwitz zeta function and the polygamma
      function, an order s with |Re s| >= 10^10 or |Im s| >= 10^7;
    * for the Barnes G-function, arguments with a part of magnitude
      >= 10^50;
    * for the zeta-type functions in the complex contexts, a huge imaginary
      part (the number of terms grows with |Im s|) or a huge real part,
      and for the Riemann xi function an imaginary part so large that the
      cancellation (xi decays like exp(-pi |Im s| / 4)) exceeds the maximal
      working precision.

    The elementary functions have their own handling of large arguments.
    re[i] and im[i] are the parts of the arguments (im may be NULL).
*/
#define DECIMAL_ARB_EXP_LIMIT 1000000

int
_decimal_arb_args_infeasible(int method, decfloat_srcptr * re, decfloat_srcptr * im, slong nargs, slong prec, gr_ctx_t ctx)
{
    slong i;

    if (_decimal_method_is_elementary(method))
        return 0;

    for (i = 0; i < 2 * nargs; i++)
    {
        decfloat_srcptr x = (i < nargs) ? re[i] : (im == NULL ? NULL : im[i - nargs]);

        if (x != NULL && !DECFLOAT_IS_SPECIAL(x))
        {
            slong E = _decfloat_sci_exp_clamped(x, ctx);

            if (E > DECIMAL_ARB_EXP_LIMIT)
                return 1;

            /* tiny parts are generally harmless (and give useful results
               near poles, for example), but not for the Barnes G-function */
            if (E < -DECIMAL_ARB_EXP_LIMIT && (method == GR_METHOD_BARNES_G || method == GR_METHOD_LOG_BARNES_G))
                return 1;
        }
    }

    switch (method)
    {
        /* huge order s */
        case GR_METHOD_POLYLOG: case GR_METHOD_HURWITZ_ZETA: case GR_METHOD_POLYGAMMA:
            if (nargs >= 1 && !DECFLOAT_IS_SPECIAL(re[0]) && _decfloat_sci_exp_clamped(re[0], ctx) >= 10)
                return 1;
            /* no Riemann-Siegel type formula: the number of terms grows
               linearly with |Im s| */
            if (nargs >= 1 && im != NULL && !DECFLOAT_IS_SPECIAL(im[0]) && _decfloat_sci_exp_clamped(im[0], ctx) >= 7)
                return 1;
            break;
        /* arb gives no finite enclosure (and spends much time) beyond this */
        case GR_METHOD_BARNES_G:
            for (i = 0; i < 2 * nargs; i++)
            {
                decfloat_srcptr x = (i < nargs) ? re[i] : (im == NULL ? NULL : im[i - nargs]);
                if (x != NULL && !DECFLOAT_IS_SPECIAL(x) && _decfloat_sci_exp_clamped(x, ctx) >= 50)
                    return 1;
            }
            break;
        default:
            break;
    }

    if (im != NULL && nargs >= 1)
    {
        switch (method)
        {
            case GR_METHOD_ZETA: case GR_METHOD_HURWITZ_ZETA: case GR_METHOD_DIRICHLET_ETA:
            case GR_METHOD_DIRICHLET_BETA: case GR_METHOD_POLYLOG: case GR_METHOD_RIEMANN_XI:
                if (!DECFLOAT_IS_SPECIAL(re[0]) && _decfloat_sci_exp_clamped(re[0], ctx) >= 30)
                    return 1;
                if (!DECFLOAT_IS_SPECIAL(im[0]))
                {
                    slong E = _decfloat_sci_exp_clamped(im[0], ctx);

                    if (E >= 15)
                        return 1;

                    /* 0.35 |Im s| digits of cancellation */
                    if (method == GR_METHOD_RIEMANN_XI && 0.35 * pow(10.0, (double) E) > (double) _decimal_ziv_wp_max(prec, 0))
                        return 1;
                }
                break;
            default:
                break;
        }
    }

    return 0;
}

int
_decimal_method_is_elementary(int method)
{
    switch (method)
    {
        case GR_METHOD_EXP: case GR_METHOD_EXPM1: case GR_METHOD_EXP2: case GR_METHOD_EXP10:
        case GR_METHOD_LOG: case GR_METHOD_LOG1P: case GR_METHOD_LOG2: case GR_METHOD_LOG10:
        case GR_METHOD_POW:
        case GR_METHOD_SIN: case GR_METHOD_COS: case GR_METHOD_TAN: case GR_METHOD_COT:
        case GR_METHOD_SEC: case GR_METHOD_CSC: case GR_METHOD_SINC:
        case GR_METHOD_SIN_PI: case GR_METHOD_COS_PI: case GR_METHOD_TAN_PI: case GR_METHOD_COT_PI:
        case GR_METHOD_SEC_PI: case GR_METHOD_CSC_PI: case GR_METHOD_SINC_PI:
        case GR_METHOD_ASIN: case GR_METHOD_ACOS: case GR_METHOD_ATAN: case GR_METHOD_ACOT:
        case GR_METHOD_ASEC: case GR_METHOD_ACSC:
        case GR_METHOD_ASIN_PI: case GR_METHOD_ACOS_PI: case GR_METHOD_ATAN_PI: case GR_METHOD_ACOT_PI:
        case GR_METHOD_ASEC_PI: case GR_METHOD_ACSC_PI:
        case GR_METHOD_SINH: case GR_METHOD_COSH: case GR_METHOD_TANH: case GR_METHOD_COTH:
        case GR_METHOD_SECH: case GR_METHOD_CSCH:
        case GR_METHOD_ASINH: case GR_METHOD_ACOSH: case GR_METHOD_ATANH: case GR_METHOD_ACOTH:
        case GR_METHOD_ASECH: case GR_METHOD_ACSCH:
        case GR_METHOD_EXP_PI_I: case GR_METHOD_LOG_PI_I:
            return 1;
        default:
            return 0;
    }
}

/* ------------------------------------------------------------------------- */
/*    Special values and tiny arguments                                     */
/* ------------------------------------------------------------------------- */

/* the leading Taylor term S of f(x) for a tiny argument, and the sign and
   size of the remaining tail: S is x, 1, +/- 1/x (if exact), x - 1 (log
   near 1), 1/(x - 1) (zeta near 1), 1/2, -1/2 or sgn(x)/2 */
enum { S_X, S_ONE, S_INVX, S_NEGINVX, S_XM1, S_INVXM1, S_HALF, S_MHALF, S_SGNHALF };

/* tail sign: fixed, or +/- the sign of x */
enum { T_POS, T_NEG, T_SGN, T_NEGSGN };

typedef struct
{
    int kind;       /* S_* */
    int tsign;      /* T_* */
    int tmul;       /* |tail| < 10^(tmul * E + tadd), E = exponent of x (or h) */
    int tadd;
    int large;      /* 0: applies to tiny x (E <= -2); 1: to large x (E >= 1) */
}
tail_info;

/* exact result kinds for _decfloat_exact_value */
enum
{
    X_NONE,
    X_GAMMA,        /* Gamma(n) = (n-1)! */
    X_RGAMMA,       /* 1/Gamma(n) */
    X_LGAMMA,       /* lgamma(1) = lgamma(2) = 0 */
    X_ZETA,         /* zeta(-n) rational, zeta(0) = -1/2 */
    X_LOG2,         /* log2(2^n) = n */
    X_LOG10,        /* log10(10^n) = n */
    X_SIN_PI, X_COS_PI, X_TAN_PI, X_COT_PI, X_CSC_PI, X_SEC_PI, X_SINC_PI,
    X_ASIN_PI, X_ACOS_PI, X_ATAN_PI, X_ACOT_PI, X_ASEC_PI, X_ACSC_PI,
    X_ACOS,         /* acos(1) = 0 */
    X_ACOSH,        /* acosh(1) = 0 */
    X_ZERO_AT_ZERO, /* f(0) = 0 (functions where arb may not return an exact zero) */
    X_ONE_AT_ZERO,  /* f(0) = 1 */
    X_BARNES_G,     /* G(n) = prod_{k<n-1} k! */
    X_LOG_BARNES_G, /* log G(1) = log G(2) = log G(3) = 0 */
    X_NUM
};


/* 1/x exactly, if it is a decimal fraction */
static int
_inv_exact(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    decfloat_t one;
    int status;
    decfloat_init(one, ctx);
    decfloat_one(one, ctx);
    status = _decfloat_div(res, one, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx);
    decfloat_clear(one, ctx);
    return status == GR_SUCCESS;
}

/* returns 1 if handled */
DECIMAL_DRIVER int
_decfloat_tiny_arg(decfloat_t res, const decfloat_t x, const tail_info * ti, slong prec, int rnd, gr_ctx_t ctx)
{
    decfloat_t Sbuf, h;
    const decfloat_struct * S;
    slong E, tail_exp;
    int sgn, tail_sign, r = 0;

    if (DECFLOAT_IS_SPECIAL(x) || prec == DECIMAL_PREC_EXACT)
        return 0;

    /* |x - 1| < 1/10 requires 1/10 <= |x| < 10: gate before the exact
       subtraction, which could be expensive for a huge or tiny x */
    if ((ti->kind == S_XM1 || ti->kind == S_INVXM1))
    {
        E = _decfloat_sci_exp_clamped(x, ctx);
        if (E != 0 && E != -1)
            return 0;
    }

    decfloat_init(Sbuf, ctx);
    decfloat_init(h, ctx);

    if (ti->kind == S_XM1 || ti->kind == S_INVXM1)
    {
        /* h = x - 1, exactly */
        decfloat_one(Sbuf, ctx);
        if (_decfloat_add(h, x, Sbuf, 1, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx) != GR_SUCCESS || DECFLOAT_IS_ZERO(h))
            goto cleanup;
        E = _decfloat_sci_exp_clamped(h, ctx);
        sgn = _decfloat_sgn(h, ctx);
    }
    else
    {
        E = _decfloat_sci_exp_clamped(x, ctx);
        sgn = _decfloat_sgn(x, ctx);
    }

    /* the tail bounds are only valid for tiny arguments; the precondition
       of _decfloat_round_with_tail makes them tiny, but check |x| < 1/10
       early to avoid pointless exact inverses */
    if (ti->large ? (E < 1) : (E > -2))
        goto cleanup;

    switch (ti->kind)
    {
        case S_X:      S = x; break;
        case S_ONE:    decfloat_one(Sbuf, ctx); S = Sbuf; break;
        case S_XM1:    S = h; break;
        case S_INVX:
        case S_NEGINVX:
            if (!_inv_exact(Sbuf, x, ctx)) goto cleanup;
            if (ti->kind == S_NEGINVX) Sbuf->m.size = -Sbuf->m.size;
            S = Sbuf;
            break;
        case S_INVXM1:
            if (!_inv_exact(Sbuf, h, ctx)) goto cleanup;
            S = Sbuf;
            break;
        case S_HALF:
        case S_MHALF:
        case S_SGNHALF:
            GR_MUST_SUCCEED(decfloat_set_si_10exp_si(Sbuf, (ti->kind == S_HALF || (ti->kind == S_SGNHALF && sgn > 0)) ? 5 : -5, -1, ctx));
            S = Sbuf;
            break;
        default: goto cleanup;
    }

    switch (ti->tsign)
    {
        case T_POS: tail_sign = 1; break;
        case T_NEG: tail_sign = -1; break;
        case T_SGN: tail_sign = sgn; break;
        default: tail_sign = -sgn; break;
    }

    /* E is clamped to +/- WORD_MAX/16 and |tmul| <= 3, so no overflow */
    tail_exp = ti->tmul * E + ti->tadd;

    r = _decfloat_round_with_tail(res, S, tail_sign, tail_exp, prec, rnd, ctx);

cleanup:
    decfloat_clear(Sbuf, ctx);
    decfloat_clear(h, ctx);
    return r;
}

/* sin(pi k / 12), cos(pi k / 12) etc. as rationals when they are rational;
   returns 1 and sets q, 0 if irrational, -1 for a pole */
static int
_trig_pi_rational(fmpq_t q, slong k, int kind)
{
    /* 2 sin(pi k/12) for k = 0..23, with 9 marking an irrational value */
    static const signed char sin2[24] = { 0, 9, 1, 9, 9, 9, 2, 9, 9, 9, 1, 9, 0, 9, -1, 9, 9, 9, -2, 9, 9, 9, -1, 9 };
    slong s, c, num, den;

    k = ((k % 24) + 24) % 24;
    s = sin2[k];
    c = sin2[(k + 6) % 24];

    switch (kind)
    {
        case X_SIN_PI: if (s == 9) return 0; num = s; den = 2; break;
        case X_COS_PI: if (c == 9) return 0; num = c; den = 2; break;
        case X_TAN_PI:
            if (c == 0) return -1;
            if (k % 6 == 3) { num = (k % 12 == 3) ? 1 : -1; den = 1; break; }   /* tan(pi/4) etc. */
            if (s == 9 || c == 9) return 0;
            num = s; den = c; break;
        case X_COT_PI:
            if (s == 0) return -1;
            if (k % 6 == 3) { num = (k % 12 == 3) ? 1 : -1; den = 1; break; }
            if (s == 9 || c == 9) return 0;
            num = c; den = s; break;
        case X_CSC_PI:
            if (s == 0) return -1;
            if (s == 9) return 0;
            num = 2; den = s; break;
        case X_SEC_PI:
            if (c == 0) return -1;
            if (c == 9) return 0;
            num = 2; den = c; break;
        default:
            return 0;
    }

    if (den < 0)
    {
        num = -num;
        den = -den;
    }
    fmpq_set_si(q, num, den);
    return 1;
}

/* Residue modulo a word-size prime p of the integer M with |x| = M 10^v,
   10 not dividing M, in time linear in the number of limbs. */
#if FLINT_BITS == 64
#define HASH_PRIME UWORD(4611686018427387847)   /* 2^62 - 57 */
#else
#define HASH_PRIME UWORD(1073741789)            /* 2^30 - 35 */
#endif

static ulong
_decfloat_mantissa_mod_p(const decfloat_t x, ulong p, ulong pinv, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    ulong B = n_mod2_preinv(DECIMAL_CTX_B(ctx), p, pinv), r = 0;
    slong i, n = FLINT_ABS(x->m.size), val;

    for (i = n - 1; i >= 0; i--)
    {
        r = n_mulmod2_preinv(r, B, p, pinv);
        r = n_addmod(r, n_mod2_preinv(x->m.d[i], p, pinv), p);
    }

    /* divide out the trailing zeros of the lowest limb */
    val = _radix_valuation_digits_1(x->m.d[0], radix);
    if (val != 0)
        r = n_mulmod2_preinv(r, n_invmod(n_powmod2_ui_preinv(10, val, p, pinv), p), p, pinv);

    return r;
}

/* Whether x = 2^n for an integer n, setting n. Cheap filters on the digits
   and a modular check come before any exact conversion. */
static int
_decfloat_is_pow2(fmpz_t n, const decfloat_t x, slong E, slong v, gr_ctx_t ctx)
{
    ulong p = HASH_PRIME, pinv = n_preinvert_limb(p), h;
    slong D = decfloat_digits(x, ctx), k, klo, khi, cand = WORD_MIN;
    ulong last = decfloat_get_digit_si(x, v, ctx);
    arf_t a;
    int r = 0;

    if (_decfloat_sgn(x, ctx) <= 0 || E == DECFLOAT_EXP_CLAMP || E == -DECFLOAT_EXP_CLAMP)
        return 0;

    if (v == 0)
    {
        /* 2^k with k >= 0: an integer with D = floor(k log10(2)) + 1 digits
           ending in 1, 2, 4, 6 or 8 */
        if (last == 0 || last == 3 || last == 5 || last == 7 || last == 9)
            return 0;
        klo = (slong) ((D - 1) * 3.3219280948873623479) - 2;
        khi = (slong) (D * 3.3219280948873623479) + 2;
        h = _decfloat_mantissa_mod_p(x, p, pinv, ctx);
        for (k = FLINT_MAX(klo, 0); k <= khi; k++)
            if (n_powmod2_ui_preinv(2, k, p, pinv) == h)
                cand = k;
    }
    else if (v < 0)
    {
        /* 2^-k = 5^k 10^-k: v = -k and M = 5^k with floor(k log10(5)) + 1
           digits, ending in 5 */
        k = -v;
        if (last != 5 || FLINT_ABS((double) (D - 1) - k * 0.69897000433601880479) > 1.0)
            return 0;
        h = _decfloat_mantissa_mod_p(x, p, pinv, ctx);
        if (n_powmod2_ui_preinv(5, k, p, pinv) == h)
            cand = -k;
    }

    if (cand == WORD_MIN)
        return 0;

    /* verify exactly */
    arf_init(a);
    if (decfloat_get_arf(a, x, ARF_PREC_EXACT, ARF_RND_DOWN, ctx) == GR_SUCCESS
        && arf_bits(a) == 1)
    {
        fmpz_sub_ui(n, ARF_EXPREF(a), 1);
        r = 1;
    }
    arf_clear(a);
    return r;
}

/* returns 1 if res was set to the exact (rounded) value, 0 if not an exact
   case, -1 on error, -2 for a domain error (pole). All tests are gated by
   the sign, the scientific exponent E, the decimal valuation v (so that
   |x| = M 10^v with 10 not dividing M) and individual digits of x, and
   only convert x to an integer or rational when it is small or when it
   is a genuine candidate. */
DECIMAL_DRIVER int
_decfloat_exact_value(decfloat_t res, const decfloat_t x, int kind, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpq_t q;
    fmpz_t n;
    slong E, v, ni = 0;
    int r = 0, sgn, is_int, is_small_int;

    if (kind == X_NONE || DECFLOAT_IS_SPECIAL(x))
    {
        if (!DECFLOAT_IS_ZERO(x))
            return 0;
        if (kind == X_ZERO_AT_ZERO || kind == X_SIN_PI || kind == X_TAN_PI || kind == X_ASIN_PI || kind == X_ATAN_PI)
            return (decfloat_zero(res, ctx) == GR_SUCCESS) ? 1 : -1;
        if (kind == X_ONE_AT_ZERO || kind == X_COS_PI || kind == X_SEC_PI || kind == X_SINC_PI)
            return (decfloat_one(res, ctx) == GR_SUCCESS) ? 1 : -1;
        if (kind == X_ACOS_PI || kind == X_ACOT_PI)
        {
            fmpq_init(q);
            fmpq_set_si(q, 1, 2);
            r = (decfloat_set_round_fmpq(res, q, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;
            fmpq_clear(q);
            return r;
        }
        if (kind == X_ZETA)
        {
            fmpq_init(q);
            fmpq_set_si(q, -1, 2);
            r = (decfloat_set_round_fmpq(res, q, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;
            fmpq_clear(q);
            return r;
        }
        if (kind == X_GAMMA || kind == X_COT_PI || kind == X_CSC_PI || kind == X_ASEC_PI || kind == X_ACSC_PI)
            return -2;
        if (kind == X_LOG2 || kind == X_LOG10)
            return (decfloat_neg_inf(res, ctx) == GR_SUCCESS) ? 1 : -2;
        if (kind == X_RGAMMA)
            return (decfloat_zero(res, ctx) == GR_SUCCESS) ? 1 : -1;
        return 0;
    }

    sgn = _decfloat_sgn(x, ctx);
    E = _decfloat_sci_exp_clamped(x, ctx);
    v = _decfloat_val10_clamped(x, ctx);
    is_int = (v >= 0);
    /* |x| < 10^(DECIMAL_WORD_DIGITS - 1) fits in a signed word */
    is_small_int = is_int && E < DECIMAL_WORD_DIGITS - 1;
    if (is_small_int)
        GR_MUST_SUCCEED(decfloat_get_si(&ni, x, ctx));

    fmpq_init(q);
    fmpz_init(n);

    switch (kind)
    {
        case X_GAMMA:
        case X_RGAMMA:
        case X_BARNES_G:
            if (!is_int) break;
            if (sgn <= 0)
            {
                /* poles of Gamma, zeros of 1/Gamma and G */
                if (kind == X_GAMMA) r = -2;
                else r = (decfloat_zero(res, ctx) == GR_SUCCESS) ? 1 : -1;
                break;
            }
            if (!is_small_int || ni > 3000) break;
            if (kind == X_BARNES_G)
            {
                slong k;
                if (ni > 200) break;
                fmpz_one(fmpq_numref(q));
                fmpz_one(fmpq_denref(q));
                for (k = 2; k + 1 < ni; k++)
                {
                    fmpz_t f;
                    fmpz_init(f);
                    fmpz_fac_ui(f, k);
                    fmpz_mul(fmpq_numref(q), fmpq_numref(q), f);
                    fmpz_clear(f);
                }
            }
            else
            {
                fmpz_fac_ui(fmpq_numref(q), ni - 1);
                fmpz_one(fmpq_denref(q));
                if (kind == X_RGAMMA)
                    fmpq_inv(q, q);
            }
            r = (decfloat_set_round_fmpq(res, q, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;
            break;

        case X_LGAMMA:
            if (is_int && sgn <= 0)
                r = -2;
            else if (is_small_int && (ni == 1 || ni == 2))
                r = (decfloat_zero(res, ctx) == GR_SUCCESS) ? 1 : -1;
            break;

        case X_LOG_BARNES_G:
            if (is_int && sgn <= 0)
                r = -2;
            else if (is_small_int && ni >= 1 && ni <= 3)
                r = (decfloat_zero(res, ctx) == GR_SUCCESS) ? 1 : -1;
            break;

        case X_ZETA:
            if (!is_int) break;
            if (is_small_int && ni == 1) { r = -2; break; }
            if (sgn < 0 && decfloat_get_digit_si(x, 0, ctx) % 2 == 0)
            {
                /* trivial zeros zeta(-2k) = 0, for any size */
                r = (decfloat_zero(res, ctx) == GR_SUCCESS) ? 1 : -1;
            }
            else if (sgn < 0 && is_small_int && ni >= -2000)
            {
                /* zeta(-m) = -B_{m+1} / (m+1), m odd */
                ulong m = -ni;
                arith_bernoulli_number(q, m + 1);
                fmpz_set_ui(n, m + 1);
                fmpq_div_fmpz(q, q, n);
                fmpq_neg(q, q);
                r = (decfloat_set_round_fmpq(res, q, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;
            }
            break;

        case X_LOG2:
            if (sgn <= 0) { r = -2; break; }
            if (_decfloat_is_pow2(n, x, E, v, ctx))
                r = (decfloat_set_round_fmpz(res, n, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;
            break;

        case X_LOG10:
            if (sgn <= 0) { r = -2; break; }
            /* x = 10^v: a single digit 1 */
            if (E == v && decfloat_get_digit_si(x, v, ctx) == 1)
            {
                decfloat_get_sci_exp(n, x, ctx);
                r = (decfloat_set_round_fmpz(res, n, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;
            }
            break;

        case X_SINC_PI:
            if (is_int)
                r = (decfloat_zero(res, ctx) == GR_SUCCESS) ? 1 : -1;
            break;

        case X_SIN_PI: case X_COS_PI: case X_TAN_PI: case X_COT_PI: case X_CSC_PI: case X_SEC_PI:
            {
                /* rational values need 12x in Z, i.e. x in Z/4 (v >= -2); the
                   value depends on 12x mod 24, i.e. on N = 100 |x| mod 200,
                   given by the digits of x at 10^0, 10^-1 and 10^-2 */
                ulong N;
                slong k;
                if (v < -2) break;
                N = decfloat_get_digit_si(x, 0, ctx) * 100 + decfloat_get_digit_si(x, -1, ctx) * 10 + decfloat_get_digit_si(x, -2, ctx);
                if (N % 25 != 0) break;
                k = (3 * (N / 25)) % 24;      /* 12 |x| mod 24 */
                if (sgn < 0) k = (24 - k) % 24;
                r = _trig_pi_rational(q, k, kind);
                if (r == 1)
                    r = (decfloat_set_round_fmpq(res, q, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;
                else if (r == -1)
                    r = -2;   /* pole */
            }
            break;

        case X_ASIN_PI: case X_ACOS_PI: case X_ATAN_PI: case X_ACOT_PI: case X_ASEC_PI: case X_ACSC_PI:
            {
                /* inputs 2x in {-4,...,4} with rational outputs, and domain
                   errors of asin, acos (|x| > 1) and asec, acsc (|x| < 1) */
                slong k2;
                if (kind == X_ASIN_PI || kind == X_ACOS_PI)
                {
                    if (E >= 1 || (E == 0 && _decfloat_cmpabs_ui(x, 1, ctx) > 0)) { r = -2; break; }
                }
                else if (kind == X_ASEC_PI || kind == X_ACSC_PI)
                {
                    if (E < 0) { r = -2; break; }
                }
                if (v < -1 || E > 0) break;
                /* |x| = (10 d_0 + d_-1) / 10 with 2|x| an integer */
                k2 = decfloat_get_digit_si(x, 0, ctx) * 10 + decfloat_get_digit_si(x, -1, ctx);
                if (k2 % 5 != 0) break;
                k2 = k2 / 5;
                if (sgn < 0) k2 = -k2;
                r = 1;
                switch (kind)
                {
                    case X_ASIN_PI:
                        if (k2 == 1) fmpq_set_si(q, 1, 6); else if (k2 == -1) fmpq_set_si(q, -1, 6);
                        else if (k2 == 2) fmpq_set_si(q, 1, 2); else if (k2 == -2) fmpq_set_si(q, -1, 2);
                        else r = (FLINT_ABS(k2) > 2) ? -2 : 0;
                        break;
                    case X_ACOS_PI:
                        if (k2 == 1) fmpq_set_si(q, 1, 3); else if (k2 == -1) fmpq_set_si(q, 2, 3);
                        else if (k2 == 2) fmpq_zero(q); else if (k2 == -2) fmpq_one(q);
                        else r = (FLINT_ABS(k2) > 2) ? -2 : 0;
                        break;
                    case X_ATAN_PI:
                    case X_ACOT_PI:
                        if (k2 == 2) fmpq_set_si(q, 1, 4); else if (k2 == -2) fmpq_set_si(q, -1, 4);
                        else r = 0;
                        break;
                    case X_ASEC_PI:
                        if (k2 == 2) fmpq_zero(q); else if (k2 == -2) fmpq_one(q);
                        else if (k2 == 4) fmpq_set_si(q, 1, 3); else if (k2 == -4) fmpq_set_si(q, 2, 3);
                        else r = (FLINT_ABS(k2) < 2) ? -2 : 0;
                        break;
                    default: /* X_ACSC_PI */
                        if (k2 == 2) fmpq_set_si(q, 1, 2); else if (k2 == -2) fmpq_set_si(q, -1, 2);
                        else if (k2 == 4) fmpq_set_si(q, 1, 6); else if (k2 == -4) fmpq_set_si(q, -1, 6);
                        else r = (FLINT_ABS(k2) < 2) ? -2 : 0;
                        break;
                }
                if (r == 1)
                    r = (decfloat_set_round_fmpq(res, q, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;
            }
            break;

        case X_ACOS:
        case X_ACOSH:
            if (is_small_int && ni == 1)
                r = (decfloat_zero(res, ctx) == GR_SUCCESS) ? 1 : -1;
            break;

        default:
            break;
    }

    fmpq_clear(q);
    fmpz_clear(n);
    return r;
}

/* ------------------------------------------------------------------------- */
/*    Function table                                                         */
/* ------------------------------------------------------------------------- */

/* Whether a finite x lies outside the real domain of an elementary function
   (cheap: a sign test and a comparison of |x| with 1). Boundary points
   (poles, branch points) are left to the other code paths. */
static int
_decfloat_outside_real_domain(const decfloat_t x, int method, gr_ctx_t ctx)
{
    int s, c;

    if (!DECFLOAT_IS_FINITE(x) || DECFLOAT_IS_ZERO(x))
        return 0;

    s = _decfloat_sgn(x, ctx);

    switch (method)
    {
        case GR_METHOD_LOG:
            return s < 0;
        case GR_METHOD_LOG1P:
            return s < 0 && _decfloat_cmpabs_ui(x, 1, ctx) > 0;
        case GR_METHOD_ASIN: case GR_METHOD_ACOS: case GR_METHOD_ATANH:
        case GR_METHOD_ERFINV:
            return _decfloat_cmpabs_ui(x, 1, ctx) > 0;
        case GR_METHOD_ASEC: case GR_METHOD_ACSC: case GR_METHOD_ACOTH:
            return _decfloat_cmpabs_ui(x, 1, ctx) < 0;
        case GR_METHOD_ACOSH:
            return s < 0 || _decfloat_cmpabs_ui(x, 1, ctx) < 0;
        case GR_METHOD_ASECH:
            return s < 0 || _decfloat_cmpabs_ui(x, 1, ctx) > 0;
        case GR_METHOD_ERFCINV:
            c = _decfloat_cmp_ui(x, 2, ctx);
            return s < 0 || c > 0;
        default:
            return 0;
    }
}

DECIMAL_DRIVER int
_decfloat_unary_fn2(decfloat_t res, const decfloat_t x, int method, int exact, const tail_info * ti, const tail_info * ti2, gr_ctx_t ctx)
{
    slong prec = DECIMAL_CTX_PREC(ctx);
    int rnd = DECIMAL_CTX_RND(ctx);
    int r;
    decfloat_srcptr a[1];
    decfloat_ptr rr[1];

    if (_decfloat_outside_real_domain(x, method, ctx))
        return GR_DOMAIN;

    r = _decfloat_exact_value(res, x, exact, prec, rnd, ctx);
    if (r == 1) return GR_SUCCESS;
    if (r == -1) return GR_UNABLE;
    if (r == -2) return GR_DOMAIN;

    if (ti != NULL)
    {
        r = _decfloat_tiny_arg(res, x, ti, prec, rnd, ctx);
        if (r == 1) return GR_SUCCESS;
        if (r == -1) return GR_UNABLE;
    }

    if (ti2 != NULL)
    {
        r = _decfloat_tiny_arg(res, x, ti2, prec, rnd, ctx);
        if (r == 1) return GR_SUCCESS;
        if (r == -1) return GR_UNABLE;
    }

    a[0] = x;
    rr[0] = res;
    return _decfloat_arb_gr(rr, 1, a, 1, 0, 0, method, ctx);
}

DECIMAL_DRIVER int
_decfloat_unary_fn(decfloat_t res, const decfloat_t x, int method, int exact, const tail_info * ti, gr_ctx_t ctx)
{
    return _decfloat_unary_fn2(res, x, method, exact, ti, NULL, ctx);
}

int
_decball_unary_fn(decball_t res, const decball_t x, int method, gr_ctx_t ctx)
{
    decball_srcptr a[1];
    decball_ptr rr[1];
    a[0] = x;
    rr[0] = res;
    return _decball_arb_gr(rr, 1, a, 1, 0, 0, method, ctx);
}

/* Functions converging exponentially fast to a representable limit for
   large arguments: tanh(x) = +/- (1 - 2 e^(-2|x|) + ...), etc. */
enum { L_TANH, L_COTH, L_ERF, L_ERFC, L_EXPM1, L_ZETA };

/* returns 1 if handled */
DECIMAL_DRIVER int
_decfloat_large_arg(decfloat_t res, const decfloat_t x, int kind, slong prec, int rnd, gr_ctx_t ctx)
{
    decfloat_t S;
    double d;
    slong tail_exp;
    int sgn, tail_sign, r = 0;

    if (DECFLOAT_IS_SPECIAL(x) || prec == DECIMAL_PREC_EXACT)
        return 0;

    sgn = _decfloat_sgn(x, ctx);

    if (kind == L_EXPM1 && sgn > 0)
        return 0;
    if (kind == L_ERFC && sgn > 0)
        return 0;
    if (kind == L_ZETA && sgn < 0)
        return 0;

    /* only large arguments */
    if (_decfloat_sci_exp_clamped(x, ctx) < 1)
        return 0;

    if (decfloat_get_d(&d, x, ctx) != GR_SUCCESS)
        d = 1e300;   /* larger than any double: the tail is negligible */
    d = fabs(d);

    /* |tail| < 10^tail_exp */
    switch (kind)
    {
        case L_TANH:
        case L_COTH:
            tail_exp = (d > 1e18) ? -DECFLOAT_EXP_CLAMP / 2 : (slong) (-0.868 * d) + 1;   /* 3 e^(-2d) < 10^(1 - 0.868 d) */
            break;
        case L_ERF:
        case L_ERFC:
            tail_exp = (d > 1e9) ? -DECFLOAT_EXP_CLAMP / 2 : (slong) (-0.434 * d * d) + 1;   /* erfc(d) < e^(-d^2) < 10^(1 - 0.434 d^2) */
            break;
        case L_ZETA:
            tail_exp = (d > 1e18) ? -DECFLOAT_EXP_CLAMP / 2 : (slong) (-0.3 * d) + 1;      /* zeta(d) - 1 < 3 * 2^-d < 10^(1 - 0.3 d) */
            break;
        default:
            tail_exp = (d > 1e18) ? -DECFLOAT_EXP_CLAMP / 2 : (slong) (-0.434 * d) + 1;       /* e^(-d) < 10^(1 - 0.434 d) */
            break;
    }

    decfloat_init(S, ctx);

    switch (kind)
    {
        case L_TANH:  decfloat_set_si(S, sgn, ctx); tail_sign = -sgn; break;   /* tanh x = sgn (1 - 2e^(-2|x|)/(1+e^(-2|x|))) */
        case L_COTH:  decfloat_set_si(S, sgn, ctx); tail_sign = sgn; break;
        case L_ERF:   decfloat_set_si(S, sgn, ctx); tail_sign = -sgn; break;
        case L_ERFC:  decfloat_set_si(S, 2, ctx); tail_sign = -1; break;       /* erfc(-d) = 2 - erfc(d) */
        case L_ZETA:  decfloat_set_si(S, 1, ctx); tail_sign = 1; break;
        default:      decfloat_set_si(S, -1, ctx); tail_sign = 1; break;       /* expm1(-d) = -1 + e^(-d) */
    }

    r = _decfloat_round_with_tail(res, S, tail_sign, tail_exp, prec, rnd, ctx);

    decfloat_clear(S, ctx);
    return r;
}

/* ------------------------------------------------------------------------- */
/*    Unary functions of real arguments                                      */
/* ------------------------------------------------------------------------- */

/* special handling before the generic path */
enum { SP_NONE, SP_LOG, SP_EXP2, SP_EXP10 };

typedef struct
{
    int method;
    signed char exact;      /* X_* */
    signed char lkind;      /* L_* for a limit at large arguments, or -1 */
    signed char special;    /* SP_* */
    signed char ntails;
    tail_info ti[2];
}
real_spec;

#define TAIL(kind, tsign, tmul, tadd) {{kind, tsign, tmul, tadd, 0}}
#define TAIL_LARGEX(tsign) {{S_INVX, tsign, -3, 1, 1}}   /* f(x) = 1/x + O(1/x^3) for large x */

/* f(x) = S + t with |t| < 10^(tmul E + tadd), for tiny x (or large x for
   the last kind), see _decfloat_tiny_arg */
static const real_spec _real_specs[] =
{
/*    method                 exact            limit    special    tails */
    {GR_METHOD_EXP,          X_NONE,          -1,      SP_NONE,   1, TAIL(S_ONE, T_SGN, 1, 2)},        /* |e^x - 1| < 2|x| */
    {GR_METHOD_EXPM1,        X_NONE,          L_EXPM1, SP_NONE,   1, TAIL(S_X, T_POS, 2, 2)},
    {GR_METHOD_EXP2,         X_NONE,          -1,      SP_EXP2,   0, {{0}}},
    {GR_METHOD_EXP10,        X_NONE,          -1,      SP_EXP10,  0, {{0}}},
    {GR_METHOD_LOG,          X_NONE,          -1,      SP_LOG,    1, TAIL(S_XM1, T_NEG, 2, 2)},        /* log(1+h) - h in (-h^2, 0) */
    {GR_METHOD_LOG1P,        X_NONE,          -1,      SP_NONE,   1, TAIL(S_X, T_NEG, 2, 2)},          /* log(1+x) - x in (-x^2, 0) */
    {GR_METHOD_LOG2,         X_LOG2,          -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_LOG10,        X_LOG10,         -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_SIN,          X_NONE,          -1,      SP_NONE,   1, TAIL(S_X, T_NEGSGN, 3, 3)},       /* |sin x - x| < |x|^3 */
    {GR_METHOD_COS,          X_NONE,          -1,      SP_NONE,   1, TAIL(S_ONE, T_NEG, 2, 2)},        /* 1 - cos x < x^2 */
    {GR_METHOD_TAN,          X_NONE,          -1,      SP_NONE,   1, TAIL(S_X, T_SGN, 3, 3)},
    {GR_METHOD_COT,          X_NONE,          -1,      SP_NONE,   1, TAIL(S_INVX, T_NEGSGN, 1, 1)},    /* cot x - 1/x = -x/3 - ... */
    {GR_METHOD_SEC,          X_NONE,          -1,      SP_NONE,   1, TAIL(S_ONE, T_POS, 2, 2)},
    {GR_METHOD_CSC,          X_NONE,          -1,      SP_NONE,   1, TAIL(S_INVX, T_SGN, 1, 1)},       /* csc x - 1/x = x/6 + ... */
    {GR_METHOD_SINC,         X_ONE_AT_ZERO,   -1,      SP_NONE,   1, TAIL(S_ONE, T_NEG, 2, 2)},        /* 1 - sinc x < x^2 */
    {GR_METHOD_SIN_PI,       X_SIN_PI,        -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_COS_PI,       X_COS_PI,        -1,      SP_NONE,   1, TAIL(S_ONE, T_NEG, 2, 3)},        /* 1 - cos(pi x) < 5 x^2 */
    {GR_METHOD_TAN_PI,       X_TAN_PI,        -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_COT_PI,       X_COT_PI,        -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_SEC_PI,       X_SEC_PI,        -1,      SP_NONE,   1, TAIL(S_ONE, T_POS, 2, 3)},        /* sec(pi x) - 1 < 5 x^2 */
    {GR_METHOD_CSC_PI,       X_CSC_PI,        -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_SINC_PI,      X_SINC_PI,       -1,      SP_NONE,   1, TAIL(S_ONE, T_NEG, 2, 3)},
    {GR_METHOD_ASIN,         X_NONE,          -1,      SP_NONE,   1, TAIL(S_X, T_SGN, 3, 3)},
    {GR_METHOD_ACOS,         X_ACOS,          -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_ATAN,         X_NONE,          -1,      SP_NONE,   1, TAIL(S_X, T_NEGSGN, 3, 3)},
    {GR_METHOD_ACOT,         X_NONE,          -1,      SP_NONE,   1, TAIL_LARGEX(T_NEGSGN)},           /* acot x = 1/x - 1/(3x^3) + ... */
    {GR_METHOD_ASEC,         X_NONE,          -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_ACSC,         X_NONE,          -1,      SP_NONE,   1, TAIL_LARGEX(T_SGN)},              /* acsc x = 1/x + 1/(6x^3) + ... */
    {GR_METHOD_ASIN_PI,      X_ASIN_PI,       -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_ACOS_PI,      X_ACOS_PI,       -1,      SP_NONE,   1, TAIL(S_HALF, T_NEGSGN, 1, 1)},    /* acos(x)/pi - 1/2 = -x/pi - ... */
    {GR_METHOD_ATAN_PI,      X_ATAN_PI,       -1,      SP_NONE,   1, {{S_SGNHALF, T_NEGSGN, -1, 0, 1}}}, /* sgn(x)/2 - 1/(pi x) + ... */
    {GR_METHOD_ACOT_PI,      X_ACOT_PI,       -1,      SP_NONE,   1, TAIL(S_SGNHALF, T_NEGSGN, 1, 1)}, /* acot(x)/pi - sgn(x)/2 = -x/pi + ... */
    {GR_METHOD_ASEC_PI,      X_ASEC_PI,       -1,      SP_NONE,   1, {{S_HALF, T_NEGSGN, -1, 0, 1}}},  /* 1/2 - 1/(pi x) - ... */
    {GR_METHOD_ACSC_PI,      X_ACSC_PI,       -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_SINH,         X_NONE,          -1,      SP_NONE,   1, TAIL(S_X, T_SGN, 3, 3)},
    {GR_METHOD_COSH,         X_NONE,          -1,      SP_NONE,   1, TAIL(S_ONE, T_POS, 2, 2)},
    {GR_METHOD_TANH,         X_NONE,          L_TANH,  SP_NONE,   1, TAIL(S_X, T_NEGSGN, 3, 3)},
    {GR_METHOD_COTH,         X_NONE,          L_COTH,  SP_NONE,   1, TAIL(S_INVX, T_SGN, 1, 1)},
    {GR_METHOD_SECH,         X_NONE,          -1,      SP_NONE,   1, TAIL(S_ONE, T_NEG, 2, 2)},
    {GR_METHOD_CSCH,         X_NONE,          -1,      SP_NONE,   1, TAIL(S_INVX, T_NEGSGN, 1, 1)},    /* csch x - 1/x = -x/6 - ... */
    {GR_METHOD_ASINH,        X_NONE,          -1,      SP_NONE,   1, TAIL(S_X, T_NEGSGN, 3, 3)},
    {GR_METHOD_ACOSH,        X_ACOSH,         -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_ATANH,        X_NONE,          -1,      SP_NONE,   1, TAIL(S_X, T_SGN, 3, 3)},
    {GR_METHOD_ACOTH,        X_NONE,          -1,      SP_NONE,   1, TAIL_LARGEX(T_SGN)},              /* acoth x = 1/x + 1/(3x^3) + ... */
    {GR_METHOD_ASECH,        X_NONE,          -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_ACSCH,        X_NONE,          -1,      SP_NONE,   1, TAIL_LARGEX(T_NEGSGN)},           /* acsch x = 1/x - 1/(6x^3) + ... */
    {GR_METHOD_LAMBERTW,     X_ZERO_AT_ZERO,  -1,      SP_NONE,   1, TAIL(S_X, T_NEG, 2, 3)},          /* W(x) - x = -x^2 + ... */
    {GR_METHOD_GAMMA,        X_GAMMA,         -1,      SP_NONE,   1, TAIL(S_INVX, T_NEG, 0, 0)},       /* Gamma(x) - 1/x = -gamma + O(x) */
    {GR_METHOD_RGAMMA,       X_RGAMMA,        -1,      SP_NONE,   1, TAIL(S_X, T_POS, 2, 2)},          /* 1/Gamma(x) - x = gamma x^2 + ... */
    {GR_METHOD_LGAMMA,       X_LGAMMA,        -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_DIGAMMA,      X_NONE,          -1,      SP_NONE,   1, TAIL(S_NEGINVX, T_NEG, 0, 0)},    /* psi(x) + 1/x = -gamma + O(x) */
    {GR_METHOD_BARNES_G,     X_BARNES_G,      -1,      SP_NONE,   1, TAIL(S_X, T_POS, 2, 2)},          /* G(x) - x = (gamma + (log(2 pi) - 1)/2) x^2 + ... */
    {GR_METHOD_LOG_BARNES_G, X_LOG_BARNES_G,  -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_ZETA,         X_ZETA,          L_ZETA,  SP_NONE,   2, {{S_INVXM1, T_POS, 0, 0, 0},     /* zeta(1+h) - 1/h = gamma + O(h) */
                                                                     {S_MHALF, T_NEGSGN, 1, 1, 0}}},  /* zeta(x) + 1/2 = -log(2 pi) x / 2 + O(x^2) */
    {GR_METHOD_ERF,          X_NONE,          L_ERF,   SP_NONE,   0, {{0}}},
    {GR_METHOD_ERFC,         X_NONE,          L_ERFC,  SP_NONE,   0, {{0}}},
    {GR_METHOD_ERFI,         X_NONE,          -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_ERFINV,       X_ZERO_AT_ZERO,  -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_ERFCINV,      X_NONE,          -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_EXP_INTEGRAL_EI, X_NONE,       -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_SIN_INTEGRAL, X_NONE,          -1,      SP_NONE,   1, TAIL(S_X, T_NEGSGN, 3, 3)},       /* Si(x) - x = -x^3/18 + ... */
    {GR_METHOD_COS_INTEGRAL, X_NONE,          -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_SINH_INTEGRAL, X_NONE,         -1,      SP_NONE,   1, TAIL(S_X, T_SGN, 3, 3)},
    {GR_METHOD_COSH_INTEGRAL, X_NONE,         -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_DILOG,        X_ZERO_AT_ZERO,  -1,      SP_NONE,   1, TAIL(S_X, T_POS, 2, 2)},          /* Li_2(x) - x = x^2/4 + ... */
    {GR_METHOD_AGM1,         X_NONE,          -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_AIRY_AI,      X_NONE,          -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_AIRY_BI,      X_NONE,          -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_AIRY_AI_PRIME, X_NONE,         -1,      SP_NONE,   0, {{0}}},
    {GR_METHOD_AIRY_BI_PRIME, X_NONE,         -1,      SP_NONE,   0, {{0}}},
    {-1, 0, 0, 0, 0, {{0}}}
};

static const real_spec *
_real_spec(int method)
{
    static signed char index[GR_METHOD_TAB_SIZE];
    static volatile int initialized = 0;
    slong i;

    if (!initialized)
    {
        for (i = 0; i < GR_METHOD_TAB_SIZE; i++)
            index[i] = -1;
        for (i = 0; _real_specs[i].method != -1; i++)
            index[_real_specs[i].method] = i;
        initialized = 1;
    }

    return (method >= 0 && method < GR_METHOD_TAB_SIZE && index[method] >= 0) ? _real_specs + index[method] : NULL;
}

/* 2^x and 10^x through pow, which detects exact cases */
static int
_decfloat_pow_base(decfloat_t res, ulong b, const decfloat_t x, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    GR_MUST_SUCCEED(decfloat_set_round_ui(t, b, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx));
    status = decfloat_pow(res, t, x, ctx);
    decfloat_clear(t, ctx);
    return status;
}

/* the unary function given by method, correctly rounded */
int
_decfloat_unary_method(decfloat_t res, const decfloat_t x, int method, gr_ctx_t ctx)
{
    const real_spec * s = _real_spec(method);
    int r;

    if (s == NULL)
        return _decfloat_unary_fn(res, x, method, X_NONE, NULL, ctx);

    if (s->special == SP_LOG && DECFLOAT_IS_ZERO(x))
        return DECIMAL_CTX_ALLOW_INF(ctx) ? decfloat_neg_inf(res, ctx) : GR_DOMAIN;
    if (s->special == SP_EXP2)
        return _decfloat_pow_base(res, 2, x, ctx);
    if (s->special == SP_EXP10)
        return _decfloat_pow_base(res, 10, x, ctx);

    if (s->lkind >= 0)
    {
        r = _decfloat_large_arg(res, x, s->lkind, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
        if (r == 1) return GR_SUCCESS;
        if (r == -1) return GR_UNABLE;
    }

    return _decfloat_unary_fn2(res, x, method, s->exact,
        (s->ntails >= 1) ? &s->ti[0] : NULL, (s->ntails >= 2) ? &s->ti[1] : NULL, ctx);
}

/* public elementary functions */
#define DEF_UNARY(name, METHOD) \
int decfloat_##name(decfloat_t res, const decfloat_t x, gr_ctx_t ctx) \
{ \
    return _decfloat_unary_method(res, x, GR_METHOD_##METHOD, ctx); \
} \
int decball_##name(decball_t res, const decball_t x, gr_ctx_t ctx) \
{ \
    return _decball_unary_fn(res, x, GR_METHOD_##METHOD, ctx); \
}

DEF_UNARY(exp, EXP)
DEF_UNARY(expm1, EXPM1)
DEF_UNARY(log, LOG)
DEF_UNARY(log1p, LOG1P)
DEF_UNARY(sin, SIN)
DEF_UNARY(cos, COS)
DEF_UNARY(tan, TAN)
DEF_UNARY(asin, ASIN)
DEF_UNARY(acos, ACOS)
DEF_UNARY(atan, ATAN)
DEF_UNARY(sinh, SINH)
DEF_UNARY(cosh, COSH)
DEF_UNARY(tanh, TANH)
DEF_UNARY(asinh, ASINH)
DEF_UNARY(acosh, ACOSH)
DEF_UNARY(atanh, ATANH)

/* sin_cos etc.: correctly rounded componentwise via the single functions
   (the tiny-argument handling is per function) */
int
_decfloat_unary2_method(decfloat_t res1, decfloat_t res2, const decfloat_t x, int method1, int method2, gr_ctx_t ctx)
{
    int status = _decfloat_unary_method(res1, x, method1, ctx);
    return status | _decfloat_unary_method(res2, x, method2, ctx);
}

int
_decball_unary2_method(decball_t res1, decball_t res2, const decball_t x, int method, gr_ctx_t ctx)
{
    decball_srcptr a[1];
    decball_ptr rr[2];
    a[0] = x; rr[0] = res1; rr[1] = res2;
    return _decball_arb_gr(rr, 2, a, 1, 0, 0, method, ctx);
}

int decfloat_sin_cos(decfloat_t res1, decfloat_t res2, const decfloat_t x, gr_ctx_t ctx) { return _decfloat_unary2_method(res1, res2, x, GR_METHOD_SIN, GR_METHOD_COS, ctx); }
int decfloat_sinh_cosh(decfloat_t res1, decfloat_t res2, const decfloat_t x, gr_ctx_t ctx) { return _decfloat_unary2_method(res1, res2, x, GR_METHOD_SINH, GR_METHOD_COSH, ctx); }
int decball_sin_cos(decball_t res1, decball_t res2, const decball_t x, gr_ctx_t ctx) { return _decball_unary2_method(res1, res2, x, GR_METHOD_SIN_COS, ctx); }
int decball_sinh_cosh(decball_t res1, decball_t res2, const decball_t x, gr_ctx_t ctx) { return _decball_unary2_method(res1, res2, x, GR_METHOD_SINH_COSH, ctx); }

/* atan2(y, x): exact values on the axes */
int
decfloat_atan2(decfloat_t res, const decfloat_t y, const decfloat_t x, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_ZERO(y) && !DECFLOAT_IS_SPECIAL(x))
    {
        if (_decfloat_sgn(x, ctx) > 0)
            return decfloat_zero(res, ctx);
        return decfloat_pi(res, ctx);
    }
    if (DECFLOAT_IS_ZERO(x) && DECFLOAT_IS_ZERO(y))
        return GR_DOMAIN;
    return _decfloat_arb_args(res, y, x, NULL, NULL, 0, DECIMAL_ARGSPEC(GR_METHOD_ATAN2, 2, 0), ctx);
}

int
decball_atan2(decball_t res, const decball_t y, const decball_t x, gr_ctx_t ctx)
{
    return _decball_arb_args(res, y, x, NULL, NULL, 0, DECIMAL_ARGSPEC(GR_METHOD_ATAN2, 2, 0), ctx);
}

/* constants */
int
_decfloat_constant(decfloat_t res, int method, gr_ctx_t ctx)
{
    decfloat_ptr rr[1];
    rr[0] = res;
    return _decfloat_arb_gr(rr, 1, NULL, 0, 0, 0, method, ctx);
}

int
_decball_constant(decball_t res, int method, gr_ctx_t ctx)
{
    decball_ptr rr[1];
    rr[0] = res;
    return _decball_arb_gr(rr, 1, NULL, 0, 0, 0, method, ctx);
}

int decfloat_pi(decfloat_t res, gr_ctx_t ctx) { return _decfloat_constant(res, GR_METHOD_PI, ctx); }
int decball_pi(decball_t res, gr_ctx_t ctx) { return _decball_constant(res, GR_METHOD_PI, ctx); }

/* ------------------------------------------------------------------------- */
/*    Integer and rational arguments                                         */
/* ------------------------------------------------------------------------- */

int
_decfloat_gamma_fmpz(decfloat_t res, const fmpz_t n, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    status = decfloat_set_round_fmpz(t, n, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
    if (status == GR_SUCCESS)
        status = _decfloat_unary_method(res, t, GR_METHOD_GAMMA, ctx);
    decfloat_clear(t, ctx);
    return status;
}

int
_decball_gamma_fmpz(decball_t res, const fmpz_t n, gr_ctx_t ctx)
{
    decball_t t;
    int status;
    decball_init(t, ctx);
    status = decball_set_fmpz(t, n, ctx);
    if (status == GR_SUCCESS)
        status = _decball_unary_fn(res, t, GR_METHOD_GAMMA, ctx);
    decball_clear(t, ctx);
    return status;
}

int
_decfloat_gamma_fmpq(decfloat_t res, const fmpq_t x, gr_ctx_t ctx)
{
    decfloat_t t;
    extra_arg extra;
    decfloat_ptr rr[1];
    int status;

    /* decimal arguments go through the usual path (exact cases etc.) */
    decfloat_init(t, ctx);
    status = decfloat_set_round_fmpq(t, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
    if (status == GR_SUCCESS)
        status = _decfloat_unary_method(res, t, GR_METHOD_GAMMA, ctx);
    decfloat_clear(t, ctx);
    if (status == GR_SUCCESS)
        return status;

    if (fmpz_sgn(fmpq_numref(x)) <= 0 && fmpz_is_one(fmpq_denref(x)))
        return GR_DOMAIN;

    extra.type = 3;
    extra.q = x;
    rr[0] = res;
    return _decfloat_arb_gr_extra(rr, 1, NULL, 0, 0, 0, &extra, GR_METHOD_GAMMA_FMPQ, ctx);
}

int
_decball_gamma_fmpq(decball_t res, const fmpq_t x, gr_ctx_t ctx)
{
    extra_arg extra;
    decball_ptr rr[1];
    extra.type = 3;
    extra.q = x;
    rr[0] = res;
    return _decball_arb_gr_extra(rr, 1, NULL, 0, 0, 0, &extra, GR_METHOD_GAMMA_FMPQ, ctx);
}

int
_decfloat_fac_ui(decfloat_t res, ulong n, gr_ctx_t ctx)
{
    if (n <= 3000)
    {
        fmpz_t t;
        int status;
        fmpz_init(t);
        fmpz_fac_ui(t, n);
        status = decfloat_set_round_fmpz(res, t, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
        fmpz_clear(t);
        return status;
    }
    else
    {
        extra_arg extra;
        decfloat_ptr rr[1];
        extra.type = 1;
        extra.u = n;
        rr[0] = res;
        return _decfloat_arb_gr_extra(rr, 1, NULL, 0, 0, 0, &extra, GR_METHOD_FAC_UI, ctx);
    }
}

int
_decball_fac_ui(decball_t res, ulong n, gr_ctx_t ctx)
{
    extra_arg extra;
    decball_ptr rr[1];
    extra.type = 1;
    extra.u = n;
    rr[0] = res;
    return _decball_arb_gr_extra(rr, 1, NULL, 0, 0, 0, &extra, GR_METHOD_FAC_UI, ctx);
}

int
_decfloat_fac_fmpz(decfloat_t res, const fmpz_t n, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    if (fmpz_sgn(n) >= 0 && fmpz_cmp_ui(n, 3000) <= 0)
        return _decfloat_fac_ui(res, fmpz_get_ui(n), ctx);
    fmpz_init(t);
    fmpz_add_ui(t, n, 1);
    status = _decfloat_gamma_fmpz(res, t, ctx);
    fmpz_clear(t);
    return status;
}

int
_decball_fac_fmpz(decball_t res, const fmpz_t n, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_add_ui(t, n, 1);
    status = _decball_gamma_fmpz(res, t, ctx);
    fmpz_clear(t);
    return status;
}

int
_decfloat_rising_ui(decfloat_t res, const decfloat_t x, ulong n, gr_ctx_t ctx)
{
    extra_arg extra;
    decfloat_srcptr a[1];
    decfloat_ptr rr[1];

    if (n == 0)
    {
        DECFLOAT_CHECK_OPERAND(x, ctx);
        return decfloat_one(res, ctx);
    }
    if (n == 1)
        return decfloat_set(res, x, ctx);

    /* exact when feasible */
    /* exact when the numerator and denominator of the exact product are
       of moderate size */
    if (!DECFLOAT_IS_SPECIAL(x) && n <= 10000
        && (double) (decfloat_digits(x, ctx) + FLINT_MIN(FLINT_ABS(_decfloat_val10_clamped(x, ctx)), 1000000)) * n <= 100000.0)
    {
        fmpq_t q, t;
        ulong k;
        int status;

        fmpq_init(q);
        fmpq_init(t);
        if (decfloat_get_fmpq(q, x, ctx) == GR_SUCCESS)
        {
            fmpq_one(t);
            for (k = 0; k < n; k++)
            {
                fmpq_mul(t, t, q);
                fmpq_add_ui(q, q, 1);
            }
            status = decfloat_set_round_fmpq(res, t, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
            fmpq_clear(q);
            fmpq_clear(t);
            return status;
        }
        fmpq_clear(q);
        fmpq_clear(t);
    }

    extra.type = 1;
    extra.u = n;
    a[0] = x;
    rr[0] = res;
    return _decfloat_arb_gr_extra(rr, 1, a, 1, 0, 0, &extra, GR_METHOD_RISING_UI, ctx);
}

int
_decball_rising_ui(decball_t res, const decball_t x, ulong n, gr_ctx_t ctx)
{
    extra_arg extra;
    decball_srcptr a[1];
    decball_ptr rr[1];
    extra.type = 1;
    extra.u = n;
    a[0] = x;
    rr[0] = res;
    return _decball_arb_gr_extra(rr, 1, a, 1, 0, 0, &extra, GR_METHOD_RISING_UI, ctx);
}

int
_decfloat_lambertw_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t k, gr_ctx_t ctx)
{
    extra_arg extra;
    decfloat_srcptr a[1];
    decfloat_ptr rr[1];
    if (fmpz_is_zero(k))
        return _decfloat_unary_method(res, x, GR_METHOD_LAMBERTW, ctx);
    extra.type = 2;
    extra.z = k;
    a[0] = x;
    rr[0] = res;
    return _decfloat_arb_gr_extra(rr, 1, a, 1, 0, 0, &extra, GR_METHOD_LAMBERTW_FMPZ, ctx);
}

int
_decball_lambertw_fmpz(decball_t res, const decball_t x, const fmpz_t k, gr_ctx_t ctx)
{
    extra_arg extra;
    decball_srcptr a[1];
    decball_ptr rr[1];
    extra.type = 2;
    extra.z = k;
    a[0] = x;
    rr[0] = res;
    return _decball_arb_gr_extra(rr, 1, a, 1, 0, 0, &extra, GR_METHOD_LAMBERTW_FMPZ, ctx);
}

POP_OPTIONS
