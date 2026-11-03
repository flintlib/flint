/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "gr.h"
#include "nfloat.h"
#include "arf/mpn_small.h"
#include "impl.h"

/*
    Division, on the quotient kernels of arf (arf/mpn_div.c): without directed
    rounding the _approx kernel (the truncation or one ulp more; exact
    for divisors with a single nonzero word, where a chain of 2/1
    divisions is used), with directed rounding the rounded kernel in the
    mode toward zero or away from zero that the sign and the context's
    direction call for. The kernels strip the zero low words of the
    mantissas and dispatch on the sizes: one and two words in registers
    (one 2/1 or two 3/2 divisions, exact with the remainder), short
    divisors by the chain, longer ones by the schoolbook approximate
    division (_approx) or flint_mpn_divhigh with a guard word and the
    exact fallback (directed).

    The vector functions keep their own paths for a fixed divisor: the
    inverse of a one- or two-word divisor is precomputed, and a longer
    one is replaced by its approximate reciprocal.
*/

/* Special values; x or y is special. Division by zero gives NaN. */
static int
_nfloat_div_special(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y) || NFLOAT_IS_ZERO(y))
        return nfloat_nan(res, ctx);

    if (NFLOAT_IS_ZERO(x))
        return nfloat_zero(res, ctx);

    if (NFLOAT_IS_INF(x))
    {
        if (NFLOAT_IS_INF(y))
            return nfloat_nan(res, ctx);

        if (NFLOAT_IS_NEG_INF(x) ^ NFLOAT_SGNBIT(y))
            return nfloat_neg_inf(res, ctx);
        else
            return nfloat_pos_inf(res, ctx);
    }

    /* finite / inf */
    return nfloat_zero(res, ctx);
}

/* The rounding mode of the kernels for a result with the sign bit
   sgnbit: the magnitude rounded away from zero when the context's
   direction asks for it, truncated otherwise. */
FLINT_FORCE_INLINE int
_nfloat_kernel_rnd(int sgnbit, gr_ctx_t ctx)
{
    return nfloat_should_round_up(sgnbit, ctx) ? ARF_RND_UP : ARF_RND_DOWN;
}

/* Writes exponent and sign of a nonzero result, rounding up the (already
   written) mantissa if needed. */
FLINT_FORCE_INLINE int
_nfloat_div_finish(nfloat_ptr res, slong exp, int sgnbit, int inexact, slong n, gr_ctx_t ctx)
{
    if (inexact && nfloat_should_round_up(sgnbit, ctx))
        NFLOAT_MANT_INCREMENT(NFLOAT_D(res), n, exp);

    NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    return GR_SUCCESS;
}

/* One- and two-word results from the exact register quotients (the
   kernels of mpn_extras/fp_small.h), rounded up if the remainder and
   the context's direction ask for it. */
FLINT_FORCE_INLINE int
_nfloat_div_1_finish(nfloat_ptr res, ulong q, ulong r, slong exp, int sgnbit, gr_ctx_t ctx)
{
    if (r != 0 && nfloat_should_round_up(sgnbit, ctx))
    {
        q++;
        if (q == 0)
        {
            q = UWORD(1) << (FLINT_BITS - 1);
            exp++;
        }
    }

    NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    NFLOAT_D(res)[0] = q;
    return GR_SUCCESS;
}

FLINT_FORCE_INLINE int
_nfloat_div_2_finish(nfloat_ptr res, ulong q1, ulong q0, int inexact, slong exp, int sgnbit, gr_ctx_t ctx)
{
    if (inexact && nfloat_should_round_up(sgnbit, ctx))
    {
        add_ssaaaa(q1, q0, q1, q0, 0, 1);
        if (q1 == 0)
        {
            q1 = UWORD(1) << (FLINT_BITS - 1);
            exp++;
        }
    }

    NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    NFLOAT_D(res)[0] = q0;
    NFLOAT_D(res)[1] = q1;
    return GR_SUCCESS;
}

FLINT_FORCE_INLINE int
_nfloat_div_1(nfloat_ptr res, ulong a, ulong b, slong exp, int sgnbit, gr_ctx_t ctx)
{
    ulong q, r, s;
    q = _flint_mpn_fdiv_1_kernel(&r, &s, a, b);
    return _nfloat_div_1_finish(res, q, r, exp + s, sgnbit, ctx);
}

FLINT_FORCE_INLINE int
_nfloat_div_2(nfloat_ptr res, ulong a1, ulong a0, ulong b1, ulong b0, slong exp, int sgnbit, gr_ctx_t ctx)
{
    ulong q1, q0, s;
    int inexact;
    inexact = _flint_mpn_fdiv_2_kernel(&q1, &q0, &s, a1, a0, b1, b0);
    return _nfloat_div_2_finish(res, q1, q0, inexact, exp + s, sgnbit, ctx);
}

/* The quotient of the normalised mantissas (a, n) and (b, bn) through
   the kernels (three words and more), with exp = exp(x) - exp(y) the
   natural exponent. */
FLINT_FORCE_INLINE int
_nfloat_div_mant(nfloat_ptr res, nn_srcptr a, nn_srcptr b, slong bn, slong exp, int sgnbit, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx), fix;

    if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
        _arf_mpn_div_approx(NFLOAT_D(res), &fix, a, n, b, bn, n * FLINT_BITS);
    else
        _arf_mpn_div(NFLOAT_D(res), &fix, a, n, b, bn, n * FLINT_BITS, _nfloat_kernel_rnd(sgnbit, ctx));

    exp += fix;
    NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    return GR_SUCCESS;
}

/******************************************************************************/
/*  Scalar functions                                                          */
/******************************************************************************/

int
nfloat_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong n, exp;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
        return _nfloat_div_special(res, x, y, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    sgnbit = NFLOAT_SGNBIT(x) ^ NFLOAT_SGNBIT(y);
    exp = NFLOAT_EXP(x) - NFLOAT_EXP(y);

    if (n == 1)
        return _nfloat_div_1(res, NFLOAT_D(x)[0], NFLOAT_D(y)[0], exp, sgnbit, ctx);
    else if (n == 2)
        return _nfloat_div_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0],
                    NFLOAT_D(y)[1], NFLOAT_D(y)[0], exp, sgnbit, ctx);
    else
        return _nfloat_div_mant(res, NFLOAT_D(x), NFLOAT_D(y), n, exp, sgnbit, ctx);
}

int
nfloat_inv(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n, exp, fix;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_INF(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    sgnbit = NFLOAT_SGNBIT(x);

    /* 1 = (1/2) 2^1 */
    if (n == 1)
        return _nfloat_div_1(res, UWORD(1) << (FLINT_BITS - 1), NFLOAT_D(x)[0], 1 - NFLOAT_EXP(x), sgnbit, ctx);
    else if (n == 2)
        return _nfloat_div_2(res, UWORD(1) << (FLINT_BITS - 1), 0, NFLOAT_D(x)[1], NFLOAT_D(x)[0], 1 - NFLOAT_EXP(x), sgnbit, ctx);

    if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
        _arf_mpn_inv_approx(NFLOAT_D(res), &fix, NFLOAT_D(x), n, n * FLINT_BITS);
    else
        _arf_mpn_inv(NFLOAT_D(res), &fix, NFLOAT_D(x), n, n * FLINT_BITS, _nfloat_kernel_rnd(sgnbit, ctx));

    exp = fix - NFLOAT_EXP(x);
    NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    return GR_SUCCESS;
}

/* Division by c = d 2^(-norm) with d normalized (dinv its inverse, or 0),
   i.e. by the nfloat with mantissa d / B and exponent FLINT_BITS - norm:
   the chain of 2/1 divisions, exact (the truncation without directed
   rounding). The scratch of the chain lives in the out-of-line function
   for two words and more, so that the one-word path has a small frame. */
static int
_nfloat_div_limb_n(nfloat_ptr res, nfloat_srcptr x, ulong d, ulong dinv, slong exp, int sgnbit, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx), fix;
    ulong scratch[2 * NFLOAT_MAX_LIMBS + 2];

    if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
        _arf_mpn_div_limb_approx(NFLOAT_D(res), &fix, NFLOAT_D(x), n, d, dinv, scratch);
    else
        _arf_mpn_div_limb(NFLOAT_D(res), &fix, NFLOAT_D(x), n, d, dinv, _nfloat_kernel_rnd(sgnbit, ctx), scratch);

    exp += fix;
    NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    return GR_SUCCESS;
}

FLINT_FORCE_INLINE int
_nfloat_div_limb(nfloat_ptr res, nfloat_srcptr x, ulong d, ulong dinv, slong norm, int sgnbit, gr_ctx_t ctx)
{
    slong exp = NFLOAT_EXP(x) - (FLINT_BITS - norm);

    if (NFLOAT_CTX_NLIMBS(ctx) == 1)
    {
        ulong q, r, s;
        q = _flint_mpn_fdiv_1_kernel_preinv(&r, &s, NFLOAT_D(x)[0], d, dinv);
        return _nfloat_div_1_finish(res, q, r, exp + s, sgnbit, ctx);
    }

    return _nfloat_div_limb_n(res, x, d, dinv, exp, sgnbit, ctx);
}

static int
_nfloat_div_ui_special(nfloat_ptr res, nfloat_srcptr x, ulong c, int csgnbit, gr_ctx_t ctx)
{
    if (c == 0 || NFLOAT_IS_NAN(x))
        return nfloat_nan(res, ctx);
    if (NFLOAT_IS_ZERO(x))
        return nfloat_zero(res, ctx);
    if (NFLOAT_IS_NEG_INF(x) ^ csgnbit)
        return nfloat_neg_inf(res, ctx);
    else
        return nfloat_pos_inf(res, ctx);
}

int
nfloat_div_ui(nfloat_ptr res, nfloat_srcptr x, ulong c, gr_ctx_t ctx)
{
    slong norm;
    ulong d, dinv;

    if (NFLOAT_IS_SPECIAL(x) || c == 0)
        return _nfloat_div_ui_special(res, x, c, 0, ctx);

    norm = flint_clz(c);
    d = c << norm;
    FLINT_MPN_INVERT_LIMB(dinv, d);
    return _nfloat_div_limb(res, x, d, dinv, norm, NFLOAT_SGNBIT(x), ctx);
}

int
nfloat_div_si(nfloat_ptr res, nfloat_srcptr x, slong c, gr_ctx_t ctx)
{
    ulong uc, d, dinv;
    slong norm;
    int csgnbit = (c < 0);

    uc = FLINT_UABS(c);

    if (NFLOAT_IS_SPECIAL(x) || c == 0)
        return _nfloat_div_ui_special(res, x, uc, csgnbit, ctx);

    norm = flint_clz(uc);
    d = uc << norm;
    FLINT_MPN_INVERT_LIMB(dinv, d);
    return _nfloat_div_limb(res, x, d, dinv, norm, NFLOAT_SGNBIT(x) ^ csgnbit, ctx);
}

/******************************************************************************/
/*  Vector functions                                                          */
/******************************************************************************/

int
_nfloat_vec_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);
    slong sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    int status = GR_SUCCESS;
    nn_ptr R = (nn_ptr) res;
    nn_srcptr X = (nn_srcptr) x;
    nn_srcptr Y = (nn_srcptr) y;

    for (i = 0; i < len; i++)
    {
        nn_ptr r = R + i * sz;
        nn_srcptr u = X + i * sz;
        nn_srcptr v = Y + i * sz;

        if (FLINT_UNLIKELY(NFLOAT_IS_SPECIAL(u) || NFLOAT_IS_SPECIAL(v)))
            status |= _nfloat_div_special(r, u, v, ctx);
        else
            status |= _nfloat_div_mant(r, NFLOAT_D(u), NFLOAT_D(v), n,
                NFLOAT_EXP(u) - NFLOAT_EXP(v), NFLOAT_SGNBIT(u) ^ NFLOAT_SGNBIT(v), ctx);
    }

    return status;
}

/* Division by a single-limb mantissa d (normalized) with exponent
   shift FLINT_BITS - norm and sign csgnbit. */
static int
_nfloat_vec_div_limb(nfloat_ptr res, nfloat_srcptr x, slong len, ulong d, slong norm, int csgnbit, gr_ctx_t ctx)
{
    slong i, sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    int status = GR_SUCCESS;
    ulong dinv;

    FLINT_MPN_INVERT_LIMB(dinv, d);

    for (i = 0; i < len; i++)
    {
        nn_ptr r = (nn_ptr) res + i * sz;
        nn_srcptr u = (nn_srcptr) x + i * sz;

        if (FLINT_UNLIKELY(NFLOAT_IS_SPECIAL(u)))
            status |= _nfloat_div_ui_special(r, u, 1, csgnbit, ctx);
        else
            status |= _nfloat_div_limb(r, u, d, dinv, norm, NFLOAT_SGNBIT(u) ^ csgnbit, ctx);
    }

    return status;
}

/* Normalizes the quotient (q, n + 1) (top limb 0 or 1) into the mantissa of
   res and finishes. q may be clobbered. */
FLINT_FORCE_INLINE int
_nfloat_div_mpn_finish(nfloat_ptr res, nn_ptr q, int inexact, slong exp, int sgnbit, slong n, gr_ctx_t ctx)
{
    if (q[n] != 0)
    {
        inexact |= (q[0] & 1);
        mpn_rshift(q, q, n + 1, 1);
        exp++;
    }

    _nfloat_copy_limbs(NFLOAT_D(res), q, n);
    return _nfloat_div_finish(res, exp, sgnbit, inexact, n, ctx);
}

/*
    Division of (x, n) by a fixed (b, n) with n >= 3 and at least two
    nonzero limbs in b, via the approximate reciprocal R = floor(B^(2n+1) / b)
    or one more, which lies in (B^(n+1), 2 B^(n+1)), so that R = B^(n+1) + R1
    with R1 < B^(n+1). For a mantissa X, the quotient with n + 1 fraction
    limbs is q = X R_exact / B^n; we compute Q = X B + H, H the high n + 1
    limbs of (X B) R1 (a lower bound with error < 2 units). Then
    Q <= X R / B^n, Q >= X R / B^n - 2 and |X R - X R_exact| / B^n < 1, so
    q lies in (Q - 1, Q + 3). The truncation of Q to n limbs is certified
    (and inexact) when the guard bits D satisfy 2 <= D <= 2^k - 4; the
    error is otherwise below 1 ulp anyway, so only directed rounding
    needs the exact fallback.
*/
static int
_nfloat_vec_div_mpn_recip(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr c, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);
    slong sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    ulong one[NFLOAT_MAX_LIMBS];
    ulong R[NFLOAT_MAX_LIMBS + 3];
    ulong Xp[NFLOAT_MAX_LIMBS + 1];
    ulong Q[NFLOAT_MAX_LIMBS + 2];
    ulong H[NFLOAT_MAX_LIMBS + 1];
    nn_srcptr b = NFLOAT_D(c);
    slong cexp = NFLOAT_EXP(c);
    int csgnbit = NFLOAT_SGNBIT(c);
    int directed = NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx);
    int status = GR_SUCCESS;

    /* (R, n + 3) = floor(B^(n-1) B^(n+2) / b) or one more */
    _nfloat_zero_limbs(one, n - 1);
    one[n - 1] = 1;
    flint_mpn_divapprox_fraction(R, one, n, b, n, n + 2);
    FLINT_ASSERT(R[n + 2] == 0 && R[n + 1] == 1);

    Xp[0] = 0;

    for (i = 0; i < len; i++)
    {
        nn_ptr r = (nn_ptr) res + i * sz;
        nn_srcptr u = (nn_srcptr) x + i * sz;
        slong exp;
        int sgnbit, certified;

        if (FLINT_UNLIKELY(NFLOAT_IS_SPECIAL(u)))
        {
            status |= _nfloat_div_special(r, u, c, ctx);
            continue;
        }

        sgnbit = NFLOAT_SGNBIT(u) ^ csgnbit;
        exp = NFLOAT_EXP(u) - cexp;

        _nfloat_copy_limbs(Xp + 1, NFLOAT_D(u), n);
        flint_mpn_mulhigh_n(H, Xp, R, n + 1);
        Q[0] = H[0];
        Q[n + 1] = mpn_add_n(Q + 1, Xp + 1, H + 1, n);

        if (!directed)
        {
            status |= _nfloat_div_mpn_finish(r, Q + 1, 0, exp, sgnbit, n, ctx);
            continue;
        }

        if (Q[n + 1] != 0)
            certified = ((Q[1] & 1) || Q[0] >= 2) && (!(Q[1] & 1) || Q[0] <= UWORD_MAX - 3);
        else
            certified = (Q[0] >= 2) && (Q[0] <= UWORD_MAX - 3);

        if (FLINT_LIKELY(certified))
            status |= _nfloat_div_mpn_finish(r, Q + 1, 1, exp, sgnbit, n, ctx);
        else
            status |= _nfloat_div_mant(r, NFLOAT_D(u), b, n, exp, sgnbit, ctx);
    }

    return status;
}

int
_nfloat_vec_div_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr c, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);
    slong sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    int status = GR_SUCCESS;
    nn_srcptr b;
    slong cexp;
    int csgnbit;

    if (len <= 1 || NFLOAT_IS_SPECIAL(c))
    {
        for (i = 0; i < len; i++)
            status |= nfloat_div((nn_ptr) res + i * sz, (nn_srcptr) x + i * sz, c, ctx);
        return status;
    }

    /* the scalar may alias an entry of res */
    {
        ulong t[NFLOAT_MAX_ALLOC];
        NFLOAT_EXP(t) = NFLOAT_EXP(c);
        NFLOAT_SGNBIT(t) = NFLOAT_SGNBIT(c);
        _nfloat_copy_limbs(NFLOAT_D(t), NFLOAT_D(c), n);
        c = t;

        b = NFLOAT_D(c);
        cexp = NFLOAT_EXP(c);
        csgnbit = NFLOAT_SGNBIT(c);

        if (n == 1)
        {
            return _nfloat_vec_div_limb(res, x, len, b[0], FLINT_BITS - cexp, csgnbit, ctx);
        }
        else if (n == 2)
        {
            ulong b1 = b[1], b0 = b[0], dinv;

            if (b0 == 0)
                return _nfloat_vec_div_limb(res, x, len, b1, FLINT_BITS - cexp, csgnbit, ctx);

            dinv = flint_mpn_preinv1(b1, b0);

            for (i = 0; i < len; i++)
            {
                nn_ptr r = (nn_ptr) res + i * sz;
                nn_srcptr u = (nn_srcptr) x + i * sz;
                ulong q1, q0, s;
                int inexact;

                if (FLINT_UNLIKELY(NFLOAT_IS_SPECIAL(u)))
                {
                    status |= _nfloat_div_special(r, u, c, ctx);
                    continue;
                }

                inexact = _flint_mpn_fdiv_2_kernel_preinv(&q1, &q0, &s,
                    NFLOAT_D(u)[1], NFLOAT_D(u)[0], b1, b0, dinv);
                status |= _nfloat_div_2_finish(r, q1, q0, inexact,
                    NFLOAT_EXP(u) - cexp + s, NFLOAT_SGNBIT(u) ^ csgnbit, ctx);
            }

            return status;
        }
        else
        {
            if (flint_mpn_zero_p(b, n - 1))
                return _nfloat_vec_div_limb(res, x, len, b[n - 1], FLINT_BITS - cexp, csgnbit, ctx);

            return _nfloat_vec_div_mpn_recip(res, x, len, c, ctx);
        }
    }
}

int
_nfloat_vec_div_scalar_ui(nfloat_ptr res, nfloat_srcptr x, slong len, ulong c, gr_ctx_t ctx)
{
    slong i, sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    int status = GR_SUCCESS;

    if (c == 0)
    {
        for (i = 0; i < len; i++)
            status |= nfloat_div_ui((nn_ptr) res + i * sz, (nn_srcptr) x + i * sz, c, ctx);
        return status;
    }

    return _nfloat_vec_div_limb(res, x, len, c << flint_clz(c), flint_clz(c), 0, ctx);
}

int
_nfloat_vec_div_scalar_si(nfloat_ptr res, nfloat_srcptr x, slong len, slong c, gr_ctx_t ctx)
{
    slong i, sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    int status = GR_SUCCESS;
    ulong uc;

    if (c == 0)
    {
        for (i = 0; i < len; i++)
            status |= nfloat_div_si((nn_ptr) res + i * sz, (nn_srcptr) x + i * sz, c, ctx);
        return status;
    }

    uc = FLINT_UABS(c);
    return _nfloat_vec_div_limb(res, x, len, uc << flint_clz(uc), flint_clz(uc), c < 0, ctx);
}
