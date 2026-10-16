/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Conversion of algebraic numbers.

    The real and imaginary parts of an algebraic number are converted
    separately. We first try a purely numerical conversion from a refined
    enclosure. Only when this cannot decide the result (Ziv's strategy
    does not terminate when the value is exactly a rounding boundary, and
    a ball should be exact when the value is exactly representable) and
    the enclosure contains a short decimal number do we check exactly
    whether the part equals that decimal number. A real irrational
    number is never a decimal number, so no exact check is needed for the
    real part of a real algebraic number of degree > 1.
*/

#include "qqbar.h"
#include "fmpq.h"
#include "arb.h"
#include "acb.h"
#include "decimal.h"
#include "gr.h"

/* Not performance-critical: optimize for size. */
PUSH_OPTIONS
OPTIMIZE_OSIZE

/* Exact sign of the real (part = 0) or imaginary (part = 1) part. */
static int
_qqbar_sgn_part(const qqbar_t x, int part)
{
    return part ? qqbar_sgn_im(x) : qqbar_sgn_re(x);
}

/* Whether the part can be rational without x being rational: the real
   part of a real algebraic number of degree > 1 is irrational. */
static int
_qqbar_part_maybe_rational(const qqbar_t x, int part)
{
    return part || !qqbar_is_real(x);
}

/* Exact test for Re(x) = q or Im(x) = q. */
static int
_qqbar_part_equal_fmpq(const qqbar_t x, int part, const fmpq_t q)
{
    qqbar_t t, u;
    int res;

    qqbar_init(t);

    if (part == 0)
    {
        /* Re(x - q) = 0; subtracting a rational is cheap */
        qqbar_sub_fmpq(t, x, q);
        res = (qqbar_sgn_re(t) == 0);
    }
    else
    {
        /* Im(x - q i) = 0 */
        qqbar_init(u);
        qqbar_i(u);
        qqbar_mul_fmpq(u, u, q);
        qqbar_sub(t, x, u);
        res = (qqbar_sgn_im(t) == 0);
        qqbar_clear(u);
    }

    qqbar_clear(t);
    return res;
}

/* Sets res to the part and returns 1 if the part is rational; returns 0
   otherwise. May compute the part exactly (slow for large degree). */
static int
_qqbar_part_get_fmpq(fmpq_t res, const qqbar_t x, int part)
{
    qqbar_t t;
    int ok;

    if (_qqbar_sgn_part(x, part) == 0)
    {
        fmpq_zero(res);
        return 1;
    }

    if (qqbar_degree(x) == 1)
    {
        qqbar_get_fmpq(res, x);
        return 1;
    }

    if (!_qqbar_part_maybe_rational(x, part))
        return 0;

    qqbar_init(t);
    if (part == 0)
        qqbar_re(t, x);
    else
        qqbar_im(t, x);
    ok = (qqbar_degree(t) == 1);
    if (ok)
        qqbar_get_fmpq(res, t);
    qqbar_clear(t);
    return ok;
}

/* Refines the enclosure z of x to about wp bits and returns a pointer
   to the requested part. */
static arb_srcptr
_qqbar_refine_part(acb_t z, const qqbar_t x, int part, slong wp)
{
    _qqbar_enclosure_raw(z, QQBAR_POLY(x), z, wp);
    return part ? acb_imagref(z) : acb_realref(z);
}

static void
_decimal_ball_ctx_init(gr_ctx_t bctx, slong prec, gr_ctx_t ctx)
{
    _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), prec, DECIMAL_RND_DOWN, 0);
    decimal_ctx_set_rad_prec(bctx, DECMAG_MAX_PREC);
}

/*
    If the ball Y contains the unique decimal number c of at most
    digits digits nearest to its midpoint, and c differs from the
    previously tested candidate last (when have_last is set), checks
    exactly whether the part of x equals c. Returns 1 (with c set) if it
    does and 0 otherwise.
*/
static int
_decimal_qqbar_candidate(decfloat_t c, decfloat_t last, int * have_last, const decball_t Y, slong digits, const qqbar_t x, int part, gr_ctx_t bctx)
{
    fmpq_t q;
    int equal = 0;

    if (decfloat_set_round(c, &Y->mid, digits, DECIMAL_RND_NEAR | DECIMAL_RND_NOLIMITS, bctx) != GR_SUCCESS)
        return 0;

    if (*have_last && _decfloat_cmp(c, last, bctx) == 0)
        return 0;

    if (!_decball_contains_decfloat(Y, c, bctx))
        return 0;

    fmpq_init(q);
    if (decfloat_get_fmpq(q, c, bctx) == GR_SUCCESS)
        equal = _qqbar_part_equal_fmpq(x, part, q);
    fmpq_clear(q);

    if (!equal)
    {
        GR_MUST_SUCCEED(decfloat_set(last, c, bctx));
        *have_last = 1;
    }

    return equal;
}

static int
_decfloat_set_round_qqbar_part(decfloat_t res, const qqbar_t x, int part, slong prec, int rnd, gr_ctx_t ctx)
{
    gr_ctx_t bctx;
    decball_t Y;
    decfloat_t c, last;
    acb_t z;
    fmpq_t q;
    slong wp;
    int status, have_last = 0, maybe_rational;

    if (_qqbar_sgn_part(x, part) == 0)
    {
        decfloat_zero(res, ctx);
        return GR_SUCCESS;
    }

    if (qqbar_degree(x) == 1 || prec == DECIMAL_PREC_EXACT)
    {
        fmpq_init(q);
        if (_qqbar_part_get_fmpq(q, x, part))
            status = decfloat_set_round_fmpq(res, q, prec, rnd, ctx);
        else
            status = GR_UNABLE;
        fmpq_clear(q);
        return status;
    }

    maybe_rational = _qqbar_part_maybe_rational(x, part);

    _decimal_ball_ctx_init(bctx, prec, ctx);
    decball_init(Y, bctx);
    decfloat_init(c, bctx);
    decfloat_init(last, bctx);
    acb_init(z);
    acb_set(z, QQBAR_ENCLOSURE(x));

    status = GR_UNABLE;

    for (wp = _decimal_digits_to_bits(prec) + 16; wp < WORD_MAX / 8; wp *= 2)
    {
        arb_srcptr v = _qqbar_refine_part(z, x, part, wp);
        int r;

        if (!arb_is_finite(v))
            continue;

        decimal_ctx_set_prec(bctx, prec + (slong) (wp * 0.30103) + 5);

        if (decball_set_arb(Y, v, bctx) != GR_SUCCESS)
            break;

        r = _decfloat_round_ball(res, Y, prec, rnd, bctx, ctx);

        if (r == -1)
            break;

        if (r == 1)
        {
            status = _decfloat_finalize(res, ctx);
            break;
        }

        /* Both endpoints may fail to round because of the exponent
           limits, which is decided by the rounding without limits. */
        if (_decfloat_round_ball(c, Y, prec, rnd | DECIMAL_RND_NOLIMITS, bctx, bctx) == 1)
        {
            status = decfloat_set_round(c, c, prec, rnd, ctx);
            if (status != GR_SUCCESS)
                break;
            status = GR_UNABLE;
        }

        /* The part is (within the enclosure) exactly on a rounding
           boundary, which is a decimal number with at most prec + 1
           digits, once the enclosure is narrow compared to the ulp. */
        if (maybe_rational && arb_rel_accuracy_bits(v) >= _decimal_digits_to_bits(prec + 1) &&
            _decimal_qqbar_candidate(c, last, &have_last, Y, prec + 1, x, part, bctx))
        {
            status = decfloat_set_round(res, c, prec, rnd, ctx);
            break;
        }
    }

    acb_clear(z);
    decfloat_clear(c, bctx);
    decfloat_clear(last, bctx);
    decball_clear(Y, bctx);
    gr_ctx_clear(bctx);
    return status;
}

static int
_decball_set_qqbar_part(decball_t res, const qqbar_t x, int part, gr_ctx_t ctx)
{
    gr_ctx_t bctx;
    decball_t Y;
    decfloat_t c, last;
    acb_t z;
    fmpq_t q;
    slong wp, prec = DECIMAL_CTX_PREC(ctx);
    int status, have_last = 0, maybe_rational;

    if (_qqbar_sgn_part(x, part) == 0)
        return decball_zero(res, ctx);

    if (qqbar_degree(x) == 1)
    {
        fmpq_init(q);
        qqbar_get_fmpq(q, x);
        status = decball_set_fmpq(res, q, ctx);
        fmpq_clear(q);
        return status;
    }

    maybe_rational = _qqbar_part_maybe_rational(x, part);

    _decimal_ball_ctx_init(bctx, prec, ctx);
    decball_init(Y, bctx);
    decfloat_init(c, bctx);
    decfloat_init(last, bctx);
    acb_init(z);
    acb_set(z, QQBAR_ENCLOSURE(x));

    status = GR_UNABLE;

    for (wp = _decimal_digits_to_bits(prec) + 16; wp < WORD_MAX / 8; wp *= 2)
    {
        arb_srcptr v = _qqbar_refine_part(z, x, part, wp);

        /* the part is nonzero, so this terminates */
        if (!arb_is_finite(v) || arb_rel_accuracy_bits(v) < _decimal_digits_to_bits(prec))
            continue;

        /* An exactly representable part has at most prec digits and is
           then the value nearest to the midpoint. */
        if (maybe_rational)
        {
            decimal_ctx_set_prec(bctx, prec + (slong) (wp * 0.30103) + 5);

            if (decball_set_arb(Y, v, bctx) == GR_SUCCESS &&
                _decimal_qqbar_candidate(c, last, &have_last, Y, prec, x, part, bctx))
            {
                status = decball_set_decfloat(res, c, ctx);
                break;
            }
        }

        status = decball_set_arb(res, v, ctx);
        break;
    }

    acb_clear(z);
    decfloat_clear(c, bctx);
    decfloat_clear(last, bctx);
    decball_clear(Y, bctx);
    gr_ctx_clear(bctx);
    return status;
}

int
decfloat_set_round_qqbar(decfloat_t res, const qqbar_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    if (qqbar_sgn_im(x) != 0)
        return GR_DOMAIN;

    return _decfloat_set_round_qqbar_part(res, x, 0, prec, rnd, ctx);
}

int
decfloat_set_qqbar(decfloat_t res, const qqbar_t x, gr_ctx_t ctx)
{
    return decfloat_set_round_qqbar(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decball_set_qqbar(decball_t res, const qqbar_t x, gr_ctx_t ctx)
{
    if (qqbar_sgn_im(x) != 0)
        return GR_DOMAIN;

    return _decball_set_qqbar_part(res, x, 0, ctx);
}

int
deccfloat_set_qqbar(deccfloat_t res, const qqbar_t x, gr_ctx_t ctx)
{
    int status;
    status = _decfloat_set_round_qqbar_part(&res->re, x, 0, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
    status |= _decfloat_set_round_qqbar_part(&res->im, x, 1, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND_IM(ctx), ctx);
    return status;
}

int
deccball_set_qqbar(deccball_t res, const qqbar_t x, gr_ctx_t ctx)
{
    int status;
    status = _decball_set_qqbar_part(&res->re, x, 0, ctx);
    status |= _decball_set_qqbar_part(&res->im, x, 1, ctx);
    return status;
}

POP_OPTIONS
