/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_ec.h"

#define AX(P) GR_EC_AFF_POINT_X(P, ctx)
#define AY(P) GR_EC_AFF_POINT_Y(P, ctx)

#define T(i) GR_ENTRY(t, i, sz)

/*
    The affine group law divides, so it is only defined over a field. A ring
    that does not advertise itself as a field is still accepted: the division
    itself reports GR_DOMAIN when it meets a non-unit.
*/
static int
_aff_require_field(gr_ec_ctx_t ctx)
{
    return (gr_ec_ctx_is_over_field(ctx) == T_FALSE) ? GR_DOMAIN : GR_SUCCESS;
}

int
gr_ec_aff_point_neg(gr_ec_aff_point_t res, const gr_ec_aff_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    if (P->is_infinity == T_TRUE)
        return gr_ec_aff_point_zero(res, ctx);
    if (P->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    if (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS)
    {
        status |= gr_set(AX(res), AX(P), R);
        status |= gr_neg(AY(res), AY(P), R);
    }
    else
    {
        gr_ptr u;

        GR_TMP_INIT(u, R);

        /* y3 = -y - a1 x - a3 */
        status |= gr_neg(u, AY(P), R);
        status |= gr_submul(u, GR_EC_A1(ctx), AX(P), R);
        status |= gr_sub(u, u, GR_EC_A3(ctx), R);
        status |= gr_set(AX(res), AX(P), R);
        status |= gr_set(AY(res), u, R);

        GR_TMP_CLEAR(u, R);
    }

    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    return status;
}

/*
    Tangent at P on y^2 + a1 x y + a3 y = x^3 + a2 x^2 + a4 x + a6:

        lambda = (3x^2 + 2 a2 x + a4 - a1 y) / (2y + a1 x + a3),
        nu     = y - lambda x,
        x3     = lambda^2 + a1 lambda - a2 - 2x,
        y3     = -(lambda + a1) x3 - nu - a3.

    Assumes P is not the point at infinity.
*/
int
_gr_ec_aff_point_dbl_long_weierstrass(gr_ec_aff_point_t res,
        const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t zero;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, 4, R);

    /* T0 = 2y + a1 x + a3 */
    status |= gr_add(T(0), AY(P), AY(P), R);
    status |= gr_addmul(T(0), GR_EC_A1(ctx), AX(P), R);
    status |= gr_add(T(0), T(0), GR_EC_A3(ctx), R);

    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR_VEC(t, 4, R);
        res->is_infinity = T_UNKNOWN;
        return status;
    }

    zero = gr_is_zero(T(0), R);

    if (zero != T_FALSE)
    {
        GR_TMP_CLEAR_VEC(t, 4, R);

        if (zero == T_TRUE)
            return gr_ec_aff_point_zero(res, ctx);

        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* T1 = lambda */
    status |= gr_sqr(T(2), AX(P), R);
    status |= gr_add(T(1), T(2), T(2), R);
    status |= gr_add(T(1), T(1), T(2), R);
    status |= gr_add(T(2), AX(P), AX(P), R);
    status |= gr_addmul(T(1), T(2), GR_EC_A2(ctx), R);
    status |= gr_add(T(1), T(1), GR_EC_A4(ctx), R);
    status |= gr_submul(T(1), GR_EC_A1(ctx), AY(P), R);

    if (status == GR_SUCCESS)
        status = gr_div(T(1), T(1), T(0), R);

    if (status == GR_SUCCESS)
    {
        /* T2 = nu = y - lambda x */
        status |= gr_mul(T(2), T(1), AX(P), R);
        status |= gr_sub(T(2), AY(P), T(2), R);

        /* T3 = x3 = lambda^2 + a1 lambda - a2 - 2x */
        status |= gr_sqr(T(3), T(1), R);
        status |= gr_addmul(T(3), GR_EC_A1(ctx), T(1), R);
        status |= gr_sub(T(3), T(3), GR_EC_A2(ctx), R);
        status |= gr_sub(T(3), T(3), AX(P), R);
        status |= gr_sub(T(3), T(3), AX(P), R);

        /* T0 = y3 = -(lambda + a1) x3 - nu - a3 */
        status |= gr_add(T(0), T(1), GR_EC_A1(ctx), R);
        status |= gr_mul(T(0), T(0), T(3), R);
        status |= gr_neg(T(0), T(0), R);
        status |= gr_sub(T(0), T(0), T(2), R);
        status |= gr_sub(T(0), T(0), GR_EC_A3(ctx), R);

        status |= gr_set(AX(res), T(3), R);
        status |= gr_set(AY(res), T(0), R);
    }

    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    GR_TMP_CLEAR_VEC(t, 4, R);

    return status;
}

/* assumes P finite and a1 = a2 = a3 = 0 */
int
_gr_ec_aff_point_dbl_short_weierstrass(gr_ec_aff_point_t res,
        const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t zero;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, 3, R);

    /* T0 = 2 y */
    status |= gr_add(T(0), AY(P), AY(P), R);

    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR_VEC(t, 3, R);
        res->is_infinity = T_UNKNOWN;
        return status;
    }

    zero = gr_is_zero(T(0), R);

    if (zero != T_FALSE)
    {
        GR_TMP_CLEAR_VEC(t, 3, R);

        if (zero == T_TRUE)
            return gr_ec_aff_point_zero(res, ctx);

        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* T1 = lambda = (3 x^2 + a4) / (2 y) */
    status |= gr_sqr(T(1), AX(P), R);
    status |= gr_add(T(2), T(1), T(1), R);
    status |= gr_add(T(1), T(1), T(2), R);
    status |= gr_add(T(1), T(1), GR_EC_A4(ctx), R);

    if (status == GR_SUCCESS)
        status = gr_div(T(1), T(1), T(0), R);

    if (status == GR_SUCCESS)
    {
        /* x3 = lambda^2 - 2 x, y3 = lambda (x - x3) - y */
        status |= gr_sqr(T(2), T(1), R);
        status |= gr_sub(T(2), T(2), AX(P), R);
        status |= gr_sub(T(2), T(2), AX(P), R);

        status |= gr_sub(T(0), AX(P), T(2), R);
        status |= gr_mul(T(0), T(0), T(1), R);
        status |= gr_sub(T(0), T(0), AY(P), R);

        status |= gr_set(AX(res), T(2), R);
        status |= gr_set(AY(res), T(0), R);
    }

    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    GR_TMP_CLEAR_VEC(t, 3, R);

    return status;
}

int
gr_ec_aff_point_dbl(gr_ec_aff_point_t res, const gr_ec_aff_point_t P,
        gr_ec_ctx_t ctx)
{
    int status = _aff_require_field(ctx);

    if (status != GR_SUCCESS)
        return status;

    if (P->is_infinity == T_TRUE)
        return gr_ec_aff_point_zero(res, ctx);
    if (P->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    if (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS)
        return _gr_ec_aff_point_dbl_short_weierstrass(res, P, ctx);

    return _gr_ec_aff_point_dbl_long_weierstrass(res, P, ctx);
}

/*
    Chord through P and Q, with lambda = (y2 - y1) / (x2 - x1) and the
    same x3, y3 as in _gr_ec_aff_point_dbl_long_weierstrass. Assumes P and
    Q are not the point at infinity.
*/
int
_gr_ec_aff_point_add_long_weierstrass(gr_ec_aff_point_t res,
        const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t zero;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, 4, R);

    /* T0 = x2 - x1 */
    status |= gr_sub(T(0), AX(Q), AX(P), R);

    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR_VEC(t, 4, R);
        res->is_infinity = T_UNKNOWN;
        return status;
    }

    zero = gr_is_zero(T(0), R);

    if (zero != T_FALSE)
    {
        truth_t same;

        if (zero == T_UNKNOWN)
        {
            GR_TMP_CLEAR_VEC(t, 4, R);
            res->is_infinity = T_UNKNOWN;
            return GR_UNABLE;
        }

        same = gr_equal(AY(P), AY(Q), R);
        GR_TMP_CLEAR_VEC(t, 4, R);

        if (same == T_TRUE)
            return _gr_ec_aff_point_dbl_long_weierstrass(res, P, ctx);
        if (same == T_FALSE)
            return gr_ec_aff_point_zero(res, ctx);

        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* T1 = lambda = (y2 - y1) / (x2 - x1) */
    status |= gr_sub(T(1), AY(Q), AY(P), R);

    if (status == GR_SUCCESS)
        status = gr_div(T(1), T(1), T(0), R);

    if (status == GR_SUCCESS)
    {
        /* T2 = nu = y1 - lambda x1 */
        status |= gr_mul(T(2), T(1), AX(P), R);
        status |= gr_sub(T(2), AY(P), T(2), R);

        /* T3 = x3 = lambda^2 + a1 lambda - a2 - x1 - x2 */
        status |= gr_sqr(T(3), T(1), R);
        status |= gr_addmul(T(3), GR_EC_A1(ctx), T(1), R);
        status |= gr_sub(T(3), T(3), GR_EC_A2(ctx), R);
        status |= gr_sub(T(3), T(3), AX(P), R);
        status |= gr_sub(T(3), T(3), AX(Q), R);

        /* T0 = y3 = -(lambda + a1) x3 - nu - a3 */
        status |= gr_add(T(0), T(1), GR_EC_A1(ctx), R);
        status |= gr_mul(T(0), T(0), T(3), R);
        status |= gr_neg(T(0), T(0), R);
        status |= gr_sub(T(0), T(0), T(2), R);
        status |= gr_sub(T(0), T(0), GR_EC_A3(ctx), R);

        status |= gr_set(AX(res), T(3), R);
        status |= gr_set(AY(res), T(0), R);
    }

    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    GR_TMP_CLEAR_VEC(t, 4, R);

    return status;
}

/* assumes P, Q finite and a1 = a2 = a3 = 0 */
int
_gr_ec_aff_point_add_short_weierstrass(gr_ec_aff_point_t res,
        const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t zero;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, 3, R);

    /* T0 = x2 - x1 */
    status |= gr_sub(T(0), AX(Q), AX(P), R);

    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR_VEC(t, 3, R);
        res->is_infinity = T_UNKNOWN;
        return status;
    }

    zero = gr_is_zero(T(0), R);

    if (zero != T_FALSE)
    {
        truth_t same;

        if (zero == T_UNKNOWN)
        {
            GR_TMP_CLEAR_VEC(t, 3, R);
            res->is_infinity = T_UNKNOWN;
            return GR_UNABLE;
        }

        same = gr_equal(AY(P), AY(Q), R);
        GR_TMP_CLEAR_VEC(t, 3, R);

        if (same == T_TRUE)
            return _gr_ec_aff_point_dbl_short_weierstrass(res, P, ctx);
        if (same == T_FALSE)
            return gr_ec_aff_point_zero(res, ctx);

        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* T1 = lambda = (y2 - y1) / (x2 - x1) */
    status |= gr_sub(T(1), AY(Q), AY(P), R);

    if (status == GR_SUCCESS)
        status = gr_div(T(1), T(1), T(0), R);

    if (status == GR_SUCCESS)
    {
        /* x3 = lambda^2 - x1 - x2, y3 = lambda (x1 - x3) - y1 */
        status |= gr_sqr(T(2), T(1), R);
        status |= gr_sub(T(2), T(2), AX(P), R);
        status |= gr_sub(T(2), T(2), AX(Q), R);

        status |= gr_sub(T(0), AX(P), T(2), R);
        status |= gr_mul(T(0), T(0), T(1), R);
        status |= gr_sub(T(0), T(0), AY(P), R);

        status |= gr_set(AX(res), T(2), R);
        status |= gr_set(AY(res), T(0), R);
    }

    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    GR_TMP_CLEAR_VEC(t, 3, R);

    return status;
}

int
gr_ec_aff_point_add(gr_ec_aff_point_t res, const gr_ec_aff_point_t P,
        const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
{
    int status = _aff_require_field(ctx);

    if (status != GR_SUCCESS)
        return status;

    if (P->is_infinity == T_TRUE)
        return gr_ec_aff_point_set(res, Q, ctx);
    if (Q->is_infinity == T_TRUE)
        return gr_ec_aff_point_set(res, P, ctx);
    if (P->is_infinity == T_UNKNOWN || Q->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    if (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS)
        return _gr_ec_aff_point_add_short_weierstrass(res, P, Q, ctx);

    return _gr_ec_aff_point_add_long_weierstrass(res, P, Q, ctx);
}

int
gr_ec_aff_point_sub(gr_ec_aff_point_t res, const gr_ec_aff_point_t P,
        const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
{
    gr_ec_aff_point_t T;
    int status;

    gr_ec_aff_point_init(T, ctx);

    status = gr_ec_aff_point_neg(T, Q, ctx);

    if (status == GR_SUCCESS)
        status = gr_ec_aff_point_add(res, P, T, ctx);

    gr_ec_aff_point_clear(T, ctx);

    return status;
}

int
_gr_ec_aff_point_mul_fmpz_binary(gr_ec_aff_point_t res,
        const gr_ec_aff_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)
{
    gr_ec_aff_point_t A, S;
    fmpz_t k;
    slong i, bits;
    int status = GR_SUCCESS;

    if (fmpz_is_zero(n))
        return gr_ec_aff_point_zero(res, ctx);

    gr_ec_aff_point_init(A, ctx);
    gr_ec_aff_point_init(S, ctx);
    fmpz_init(k);
    fmpz_abs(k, n);

    if (fmpz_sgn(n) < 0)
        status |= gr_ec_aff_point_neg(S, P, ctx);
    else
        status |= gr_ec_aff_point_set(S, P, ctx);

    bits = fmpz_bits(k);

    for (i = bits - 1; i >= 0 && status == GR_SUCCESS; i--)
    {
        status |= gr_ec_aff_point_dbl(A, A, ctx);

        if (fmpz_tstbit(k, i))
            status |= gr_ec_aff_point_add(A, A, S, ctx);
    }

    if (status == GR_SUCCESS)
        status = gr_ec_aff_point_set(res, A, ctx);
    else
        res->is_infinity = T_UNKNOWN;

    fmpz_clear(k);
    gr_ec_aff_point_clear(A, ctx);
    gr_ec_aff_point_clear(S, ctx);

    return status;
}

/*
    Multiply in Jacobian coordinates, where the ladder needs no inversion,
    and convert back at the end. The conversion in is free, since an affine
    point is a Jacobian point with Z = 1.
*/
int
gr_ec_aff_point_mul_fmpz(gr_ec_aff_point_t res, const gr_ec_aff_point_t P,
        const fmpz_t n, gr_ec_ctx_t ctx)
{
    gr_ec_jac_point_t J;
    int status = _aff_require_field(ctx);

    if (status != GR_SUCCESS)
        return status;

    if (fmpz_is_zero(n) || P->is_infinity == T_TRUE)
        return gr_ec_aff_point_zero(res, ctx);

    if (P->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    gr_ec_jac_point_init(J, ctx);

    status = gr_ec_jac_point_set_aff_point(J, P, ctx);

    if (status == GR_SUCCESS)
        status = gr_ec_jac_point_mul_fmpz(J, J, n, ctx);

    if (status == GR_SUCCESS)
        status = gr_ec_aff_point_set_jac_point(res, J, ctx);
    else
        res->is_infinity = T_UNKNOWN;

    gr_ec_jac_point_clear(J, ctx);

    return status;
}

int
gr_ec_aff_point_mul_ui(gr_ec_aff_point_t res, const gr_ec_aff_point_t P,
        ulong n, gr_ec_ctx_t ctx)
{
    fmpz_t t;
    int status;

    fmpz_init_set_ui(t, n);
    status = gr_ec_aff_point_mul_fmpz(res, P, t, ctx);
    fmpz_clear(t);

    return status;
}

int
gr_ec_aff_point_mul_si(gr_ec_aff_point_t res, const gr_ec_aff_point_t P,
        slong n, gr_ec_ctx_t ctx)
{
    fmpz_t t;
    int status;

    fmpz_init(t);
    fmpz_set_si(t, n);
    status = gr_ec_aff_point_mul_fmpz(res, P, t, ctx);
    fmpz_clear(t);

    return status;
}

/*
    2^k P, by repeated doubling. A negative k would be point halving:
    that is a real operation on a curve, but it is not single valued (a
    point has up to four halves) and computing it needs the division
    polynomials, so it is out of the domain here.
*/
int
gr_ec_aff_point_mul_2exp_si(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, slong k,
        gr_ec_ctx_t ctx)
{
    int status;
    slong i;

    if (k < 0)
        return GR_DOMAIN;

    if (_aff_require_field(ctx) != GR_SUCCESS)
        return GR_DOMAIN;

    status = gr_ec_aff_point_set(res, P, ctx);

    for (i = 0; i < k && status == GR_SUCCESS; i++)
        status = gr_ec_aff_point_dbl(res, res, ctx);

    return status;
}

int
gr_ec_aff_point_mul_2exp_fmpz(gr_ec_aff_point_t res, const gr_ec_aff_point_t P,
        const fmpz_t k, gr_ec_ctx_t ctx)
{
    if (fmpz_sgn(k) < 0)
        return GR_DOMAIN;

    if (_aff_require_field(ctx) != GR_SUCCESS)
        return GR_DOMAIN;

    /* 2^k O = O for every k, including one too large to iterate over */
    if (gr_ec_aff_point_is_inf(P, ctx) == T_TRUE)
        return gr_ec_aff_point_zero(res, ctx);

    if (!fmpz_fits_si(k))
        return GR_UNABLE;

    return gr_ec_aff_point_mul_2exp_si(res, P, fmpz_get_si(k), ctx);
}
