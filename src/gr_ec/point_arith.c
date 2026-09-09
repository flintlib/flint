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

#define PX(P) GR_EC_POINT_X(P, ctx)
#define PY(P) GR_EC_POINT_Y(P, ctx)
#define PZ(P) GR_EC_POINT_Z(P, ctx)

#define T(i) GR_ENTRY(t, i, sz)

int
gr_ec_point_neg(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, R);

    /* Y3 = -Y - a1 X - a3 Z */
    status |= gr_neg(t, PY(P), R);
    status |= gr_submul(t, GR_EC_A1(ctx), PX(P), R);
    status |= gr_submul(t, GR_EC_A3(ctx), PZ(P), R);

    status |= gr_set(PX(res), PX(P), R);
    status |= gr_set(PZ(res), PZ(P), R);
    status |= gr_set(PY(res), t, R);

    GR_TMP_CLEAR(t, R);

    return status;
}

/*
    Doubling in the general Weierstrass form. In affine coordinates the
    tangent slope is lambda = u/v with

        u = 3x^2 + 2 a2 x + a4 - a1 y,  v = 2y + a1 x + a3,

    and (x3, y3) = (lambda^2 + a1 lambda - a2 - 2x, -(lambda + a1) x3 - nu - a3)
    with nu = y - lambda x. Homogenizing with u = M, v = Z S and clearing
    denominators gives the formulas below; they reduce to dbl-1998-cmo-2
    when a1 = a2 = a3 = 0. Assumes P is not the point at infinity.
*/
int
_gr_ec_point_dbl_long_weierstrass(gr_ec_point_t res, const gr_ec_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t zero;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, 9, R);

    /* T0 = S = 2 Y + a1 X + a3 Z */
    status |= gr_add(T(0), PY(P), PY(P), R);
    status |= gr_addmul(T(0), GR_EC_A1(ctx), PX(P), R);
    status |= gr_addmul(T(0), GR_EC_A3(ctx), PZ(P), R);

    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR_VEC(t, 9, R);
        return status;
    }

    zero = gr_is_zero(T(0), R);

    if (zero != T_FALSE)
    {
        GR_TMP_CLEAR_VEC(t, 9, R);
        return (zero == T_TRUE) ? gr_ec_point_zero(res, ctx) : GR_UNABLE;
    }

    /* T1 = u = 3 X^2 + 2 a2 X Z + a4 Z^2 - a1 Y Z, T2 = v = Z S */
    status |= gr_sqr(T(3), PX(P), R);
    status |= gr_add(T(1), T(3), T(3), R);
    status |= gr_add(T(1), T(1), T(3), R);
    status |= gr_mul(T(3), PX(P), PZ(P), R);
    status |= gr_add(T(4), T(3), T(3), R);
    status |= gr_addmul(T(1), T(4), GR_EC_A2(ctx), R);
    status |= gr_sqr(T(4), PZ(P), R);
    status |= gr_addmul(T(1), T(4), GR_EC_A4(ctx), R);
    status |= gr_mul(T(4), PY(P), PZ(P), R);
    status |= gr_submul(T(1), T(4), GR_EC_A1(ctx), R);

    status |= gr_mul(T(2), PZ(P), T(0), R);

    /* T3 = vv, T4 = vvv */
    status |= gr_sqr(T(3), T(2), R);
    status |= gr_mul(T(4), T(3), T(2), R);

    /* T5 = A = u^2 Z + a1 u v Z - a2 vv Z - 2 X vv */
    status |= gr_sqr(T(5), T(1), R);
    status |= gr_mul(T(5), T(5), PZ(P), R);
    status |= gr_mul(T(6), T(1), T(2), R);
    status |= gr_mul(T(6), T(6), PZ(P), R);
    status |= gr_addmul(T(5), T(6), GR_EC_A1(ctx), R);
    status |= gr_mul(T(6), T(3), PZ(P), R);
    status |= gr_submul(T(5), T(6), GR_EC_A2(ctx), R);
    status |= gr_mul(T(6), PX(P), T(3), R);
    status |= gr_sub(T(5), T(5), T(6), R);
    status |= gr_sub(T(5), T(5), T(6), R);

    /* T7 = Y3 = -(u + a1 v) A - (Y v - u X) vv - a3 vvv Z */
    status |= gr_mul(T(6), GR_EC_A1(ctx), T(2), R);
    status |= gr_add(T(6), T(6), T(1), R);
    status |= gr_mul(T(7), T(6), T(5), R);
    status |= gr_neg(T(7), T(7), R);

    status |= gr_mul(T(6), PY(P), T(2), R);
    status |= gr_submul(T(6), T(1), PX(P), R);
    status |= gr_mul(T(6), T(6), T(3), R);
    status |= gr_sub(T(7), T(7), T(6), R);

    /* T8 = Z3 = vvv Z */
    status |= gr_mul(T(8), T(4), PZ(P), R);
    status |= gr_submul(T(7), T(8), GR_EC_A3(ctx), R);

    /* T6 = X3 = v A */
    status |= gr_mul(T(6), T(2), T(5), R);

    status |= gr_set(PX(res), T(6), R);
    status |= gr_set(PY(res), T(7), R);
    status |= gr_set(PZ(res), T(8), R);

    GR_TMP_CLEAR_VEC(t, 9, R);

    return status;
}

/* dbl-1998-cmo-2; assumes P is not the point at infinity */
int
_gr_ec_point_dbl_short_weierstrass(gr_ec_point_t res, const gr_ec_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t zero;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, 10, R);

    /* T0 = XX, T2 = w = a4 ZZ + 3 XX, T3 = s = 2 Y Z */
    status |= gr_sqr(T(0), PX(P), R);
    status |= gr_sqr(T(1), PZ(P), R);
    status |= gr_mul(T(2), T(1), GR_EC_A4(ctx), R);
    status |= gr_add(T(1), T(0), T(0), R);
    status |= gr_add(T(1), T(1), T(0), R);
    status |= gr_add(T(2), T(2), T(1), R);

    status |= gr_mul(T(3), PY(P), PZ(P), R);
    status |= gr_add(T(3), T(3), T(3), R);

    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR_VEC(t, 10, R);
        return status;
    }

    zero = gr_is_zero(T(3), R);

    if (zero != T_FALSE)
    {
        GR_TMP_CLEAR_VEC(t, 10, R);
        return (zero == T_TRUE) ? gr_ec_point_zero(res, ctx) : GR_UNABLE;
    }

    /* T4 = ss, T5 = sss, T6 = Rr = Y s, T7 = RR, T8 = B, T9 = h */
    status |= gr_sqr(T(4), T(3), R);
    status |= gr_mul(T(5), T(4), T(3), R);
    status |= gr_mul(T(6), PY(P), T(3), R);
    status |= gr_sqr(T(7), T(6), R);

    /* B = (X + Rr)^2 - XX - RR */
    status |= gr_add(T(8), PX(P), T(6), R);
    status |= gr_sqr(T(8), T(8), R);
    status |= gr_sub(T(8), T(8), T(0), R);
    status |= gr_sub(T(8), T(8), T(7), R);

    /* h = w^2 - 2 B */
    status |= gr_sqr(T(9), T(2), R);
    status |= gr_sub(T(9), T(9), T(8), R);
    status |= gr_sub(T(9), T(9), T(8), R);

    /* Y3 = w (B - h) - 2 RR */
    status |= gr_sub(T(8), T(8), T(9), R);
    status |= gr_mul(T(8), T(8), T(2), R);
    status |= gr_sub(T(8), T(8), T(7), R);
    status |= gr_sub(T(8), T(8), T(7), R);

    /* X3 = h s */
    status |= gr_mul(T(9), T(9), T(3), R);

    status |= gr_set(PX(res), T(9), R);
    status |= gr_set(PY(res), T(8), R);
    status |= gr_set(PZ(res), T(5), R);

    GR_TMP_CLEAR_VEC(t, 10, R);

    return status;
}

int
gr_ec_point_dbl(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)
{
    truth_t inf = gr_ec_point_is_inf(P, ctx);

    if (inf == T_TRUE)
        return gr_ec_point_zero(res, ctx);
    if (inf == T_UNKNOWN)
        return GR_UNABLE;

    if (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS)
        return _gr_ec_point_dbl_short_weierstrass(res, P, ctx);

    return _gr_ec_point_dbl_long_weierstrass(res, P, ctx);
}

/*
    Addition in the general Weierstrass form: the chord slope is
    lambda = u/v with u = Y2 Z1 - Y1 Z2 and v = X2 Z1 - X1 Z2, and the
    affine formulas of _gr_ec_point_dbl_long_weierstrass homogenize to the
    expressions below. They reduce to add-1998-cmo-2 when a1 = a2 = a3 = 0.
    Assumes P and Q are not the point at infinity.
*/
int
_gr_ec_point_add_long_weierstrass(gr_ec_point_t res, const gr_ec_point_t P,
        const gr_ec_point_t Q, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t zero;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, 11, R);

    /* T0 = X1 Z2, T1 = X2 Z1, T2 = v, T3 = Y1 Z2, T4 = u */
    status |= gr_mul(T(0), PX(P), PZ(Q), R);
    status |= gr_mul(T(1), PX(Q), PZ(P), R);
    status |= gr_sub(T(2), T(1), T(0), R);
    status |= gr_mul(T(3), PY(P), PZ(Q), R);
    status |= gr_mul(T(4), PY(Q), PZ(P), R);
    status |= gr_sub(T(4), T(4), T(3), R);

    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR_VEC(t, 11, R);
        return status;
    }

    zero = gr_is_zero(T(2), R);

    if (zero != T_FALSE)
    {
        truth_t same;

        if (zero == T_UNKNOWN)
        {
            GR_TMP_CLEAR_VEC(t, 11, R);
            return GR_UNABLE;
        }

        same = gr_is_zero(T(4), R);
        GR_TMP_CLEAR_VEC(t, 11, R);

        if (same == T_TRUE)
            return _gr_ec_point_dbl_long_weierstrass(res, P, ctx);
        if (same == T_FALSE)
            return gr_ec_point_zero(res, ctx);

        return GR_UNABLE;
    }

    /* T5 = w = Z1 Z2, T6 = vv, T7 = vvv */
    status |= gr_mul(T(5), PZ(P), PZ(Q), R);
    status |= gr_sqr(T(6), T(2), R);
    status |= gr_mul(T(7), T(6), T(2), R);

    /* T8 = A = u^2 w + a1 u v w - a2 vv w - vv (X1 Z2 + X2 Z1) */
    status |= gr_sqr(T(8), T(4), R);
    status |= gr_mul(T(8), T(8), T(5), R);
    status |= gr_mul(T(9), T(4), T(2), R);
    status |= gr_mul(T(9), T(9), T(5), R);
    status |= gr_addmul(T(8), T(9), GR_EC_A1(ctx), R);
    status |= gr_mul(T(9), T(6), T(5), R);
    status |= gr_submul(T(8), T(9), GR_EC_A2(ctx), R);
    status |= gr_add(T(9), T(0), T(1), R);
    status |= gr_mul(T(9), T(9), T(6), R);
    status |= gr_sub(T(8), T(8), T(9), R);

    /* T10 = Y3 = -(u + a1 v) A - (Y1 v - u X1) vv Z2 - a3 vvv w */
    status |= gr_mul(T(9), GR_EC_A1(ctx), T(2), R);
    status |= gr_add(T(9), T(9), T(4), R);
    status |= gr_mul(T(10), T(9), T(8), R);
    status |= gr_neg(T(10), T(10), R);

    status |= gr_mul(T(9), PY(P), T(2), R);
    status |= gr_submul(T(9), T(4), PX(P), R);
    status |= gr_mul(T(9), T(9), T(6), R);
    status |= gr_mul(T(9), T(9), PZ(Q), R);
    status |= gr_sub(T(10), T(10), T(9), R);

    /* T7 = Z3 = vvv w */
    status |= gr_mul(T(7), T(7), T(5), R);
    status |= gr_submul(T(10), T(7), GR_EC_A3(ctx), R);

    /* T9 = X3 = v A */
    status |= gr_mul(T(9), T(2), T(8), R);

    status |= gr_set(PX(res), T(9), R);
    status |= gr_set(PY(res), T(10), R);
    status |= gr_set(PZ(res), T(7), R);

    GR_TMP_CLEAR_VEC(t, 11, R);

    return status;
}

/* add-1998-cmo-2; assumes P and Q are not the point at infinity */
int
_gr_ec_point_add_short_weierstrass(gr_ec_point_t res, const gr_ec_point_t P,
        const gr_ec_point_t Q, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t zero;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, 9, R);

    /* T0 = Y1 Z2, T1 = X1 Z2, T2 = Z1 Z2, T3 = u, T4 = v */
    status |= gr_mul(T(0), PY(P), PZ(Q), R);
    status |= gr_mul(T(1), PX(P), PZ(Q), R);
    status |= gr_mul(T(2), PZ(P), PZ(Q), R);

    status |= gr_mul(T(3), PY(Q), PZ(P), R);
    status |= gr_sub(T(3), T(3), T(0), R);

    status |= gr_mul(T(4), PX(Q), PZ(P), R);
    status |= gr_sub(T(4), T(4), T(1), R);

    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR_VEC(t, 9, R);
        return status;
    }

    zero = gr_is_zero(T(4), R);

    if (zero != T_FALSE)
    {
        truth_t same;

        if (zero == T_UNKNOWN)
        {
            GR_TMP_CLEAR_VEC(t, 9, R);
            return GR_UNABLE;
        }

        /* same x-coordinate: P = Q or P = -Q */
        same = gr_is_zero(T(3), R);
        GR_TMP_CLEAR_VEC(t, 9, R);

        if (same == T_TRUE)
            return _gr_ec_point_dbl_short_weierstrass(res, P, ctx);
        if (same == T_FALSE)
            return gr_ec_point_zero(res, ctx);

        return GR_UNABLE;
    }

    /* T5 = vv, T6 = vvv, T7 = Rr = vv X1Z2, T8 = A */
    status |= gr_sqr(T(5), T(4), R);
    status |= gr_mul(T(6), T(5), T(4), R);
    status |= gr_mul(T(7), T(5), T(1), R);

    status |= gr_sqr(T(8), T(3), R);
    status |= gr_mul(T(8), T(8), T(2), R);
    status |= gr_sub(T(8), T(8), T(6), R);
    status |= gr_sub(T(8), T(8), T(7), R);
    status |= gr_sub(T(8), T(8), T(7), R);

    /* Y3 = u (Rr - A) - vvv Y1Z2 */
    status |= gr_sub(T(7), T(7), T(8), R);
    status |= gr_mul(T(7), T(7), T(3), R);
    status |= gr_mul(T(0), T(6), T(0), R);
    status |= gr_sub(T(7), T(7), T(0), R);

    /* X3 = v A, Z3 = vvv Z1Z2 */
    status |= gr_mul(T(8), T(4), T(8), R);
    status |= gr_mul(T(6), T(6), T(2), R);

    status |= gr_set(PX(res), T(8), R);
    status |= gr_set(PY(res), T(7), R);
    status |= gr_set(PZ(res), T(6), R);

    GR_TMP_CLEAR_VEC(t, 9, R);

    return status;
}

int
gr_ec_point_add(gr_ec_point_t res, const gr_ec_point_t P,
        const gr_ec_point_t Q, gr_ec_ctx_t ctx)
{
    truth_t inf;

    inf = gr_ec_point_is_inf(P, ctx);
    if (inf == T_TRUE)
        return gr_ec_point_set(res, Q, ctx);
    if (inf == T_UNKNOWN)
        return GR_UNABLE;

    inf = gr_ec_point_is_inf(Q, ctx);
    if (inf == T_TRUE)
        return gr_ec_point_set(res, P, ctx);
    if (inf == T_UNKNOWN)
        return GR_UNABLE;

    if (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS)
        return _gr_ec_point_add_short_weierstrass(res, P, Q, ctx);

    return _gr_ec_point_add_long_weierstrass(res, P, Q, ctx);
}

int
gr_ec_point_sub(gr_ec_point_t res, const gr_ec_point_t P,
        const gr_ec_point_t Q, gr_ec_ctx_t ctx)
{
    gr_ec_point_t T;
    int status;

    gr_ec_point_init(T, ctx);

    status = gr_ec_point_neg(T, Q, ctx);

    if (status == GR_SUCCESS)
        status = gr_ec_point_add(res, P, T, ctx);

    gr_ec_point_clear(T, ctx);

    return status;
}

int
_gr_ec_point_mul_fmpz_binary(gr_ec_point_t res, const gr_ec_point_t P,
        const fmpz_t n, gr_ec_ctx_t ctx)
{
    gr_ec_point_t A, S;
    fmpz_t k;
    slong i, bits;
    int status = GR_SUCCESS;

    if (fmpz_is_zero(n))
        return gr_ec_point_zero(res, ctx);

    gr_ec_point_init(A, ctx);
    gr_ec_point_init(S, ctx);
    fmpz_init(k);
    fmpz_abs(k, n);

    if (fmpz_sgn(n) < 0)
        status |= gr_ec_point_neg(S, P, ctx);
    else
        status |= gr_ec_point_set(S, P, ctx);

    bits = fmpz_bits(k);

    for (i = bits - 1; i >= 0 && status == GR_SUCCESS; i--)
    {
        status |= gr_ec_point_dbl(A, A, ctx);

        if (fmpz_tstbit(k, i))
            status |= gr_ec_point_add(A, A, S, ctx);
    }

    if (status == GR_SUCCESS)
        status = gr_ec_point_set(res, A, ctx);

    fmpz_clear(k);
    gr_ec_point_clear(A, ctx);
    gr_ec_point_clear(S, ctx);

    return status;
}

/* Jacobian coordinates are cheaper and the conversions are inversion-free */
int
gr_ec_point_mul_fmpz(gr_ec_point_t res, const gr_ec_point_t P, const fmpz_t n,
        gr_ec_ctx_t ctx)
{
    gr_ec_jac_point_t J;
    int status;

    gr_ec_jac_point_init(J, ctx);

    status = gr_ec_jac_point_set_point(J, P, ctx);

    if (status == GR_SUCCESS)
        status = gr_ec_jac_point_mul_fmpz(J, J, n, ctx);

    if (status == GR_SUCCESS)
        status = gr_ec_point_set_jac_point(res, J, ctx);

    gr_ec_jac_point_clear(J, ctx);

    return status;
}

int
gr_ec_point_mul_ui(gr_ec_point_t res, const gr_ec_point_t P, ulong n,
        gr_ec_ctx_t ctx)
{
    fmpz_t t;
    int status;

    fmpz_init_set_ui(t, n);
    status = gr_ec_point_mul_fmpz(res, P, t, ctx);
    fmpz_clear(t);

    return status;
}

int
gr_ec_point_mul_si(gr_ec_point_t res, const gr_ec_point_t P, slong n,
        gr_ec_ctx_t ctx)
{
    fmpz_t t;
    int status;

    fmpz_init(t);
    fmpz_set_si(t, n);
    status = gr_ec_point_mul_fmpz(res, P, t, ctx);
    fmpz_clear(t);

    return status;
}
